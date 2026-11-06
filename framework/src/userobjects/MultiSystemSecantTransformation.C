//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultiSystemSecantTransformation.h"
#include "AllLocalDofIndicesThread.h"
#include "FEProblemBase.h"
#include "LinearSystem.h"

#include "libmesh/numeric_vector.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

registerMooseObject("MooseApp", MultiSystemSecantTransformation);

InputParameters
MultiSystemSecantTransformation::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Applies a secant-type acceleration (per-DOF secant, Aitken's dynamic relaxation or "
      "Anderson acceleration with depth 1) to selected variables right after the solve of their "
      "system in each multi-system fixed point iteration.");
  params.addRequiredParam<std::vector<VariableName>>(
      "variables", "Solver variables to transform. They must all belong to the same solver system.");
  params.addParam<Real>(
      "initial_relaxation_factor",
      1.0,
      "Relaxation factor w_0 of the first iteration, x + w_0 (f(x) - x), when no secant "
      "information is available yet. Only used with 'per_dof' and 'aitken'; Anderson uses "
      "'anderson_mixing' instead");
  // The mixing (damping) parameter of Anderson acceleration is taken in (0, 1], 1 being the
  // undamped update (Walker and Ni, SIAM J. Numer. Anal. 49(4), 2011; Evans et al., SIAM J.
  // Numer. Anal. 58(1), 2020)
  params.addRangeCheckedParam<Real>(
      "anderson_mixing",
      1.0,
      "anderson_mixing>0 & anderson_mixing<=1",
      "Mixing parameter beta of the Anderson acceleration, applied at every iteration "
      "including the first one (only used with method = anderson)");
  MooseEnum method("per_dof aitken anderson", "per_dof");
  params.addParam<MooseEnum>(
      "method",
      method,
      "Acceleration to use: 'per_dof' applies an independent secant update to each degree of "
      "freedom, 'aitken' computes a single dynamic relaxation factor for all degrees of freedom, "
      "'anderson' is Anderson acceleration with depth 1 (vector secant minimizing the "
      "linearized fixed point residual)");
  MooseEnum aitken_formula("recursive crossed_secant", "recursive");
  params.addParam<MooseEnum>(
      "aitken_formula",
      aitken_formula,
      "Formula of the Aitken relaxation factor (only used with method = aitken). "
      "With dx = x_k - x_{k-1} and dr = r_k - r_{k-1}: 'recursive' is the Irons-Tuck formula "
      "w_k = -w_{k-1} r_{k-1}.dr / |dr|^2, which relies on the previous factor; "
      "'crossed_secant' is its non-recursive form w_k = -dx.dr / |dr|^2, which only uses the "
      "last two iterates and is therefore independent of how x_k was obtained");
  params.addParamNamesToGroup("aitken_formula anderson_mixing", "Advanced");

  params.addParam<bool>("verbose",
                        false,
                        "Print at each iteration the coefficients of the update, and the effective "
                        "relaxation factor w_eff = (x_{k+1} - x_k).r_k / |r_k|^2");

  // The transformation is applied right after the solve of the system holding the variables, and
  // the iteration history is reset at the end of the loop
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum = {EXEC_MULTISYSTEM_FIXED_POINT_TRANSFORMATION, EXEC_MULTISYSTEM_FIXED_POINT_END};
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

MultiSystemSecantTransformation::MultiSystemSecantTransformation(
    const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _sys(_fe_problem.getVariable(0, getParam<std::vector<VariableName>>("variables")[0]).sys()),
    _var_names(getParam<std::vector<VariableName>>("variables")),
    _initial_relaxation_factor(getParam<Real>("initial_relaxation_factor")),
    _anderson_mixing(getParam<Real>("anderson_mixing")),
    _verbose(getParam<bool>("verbose")),
    _method(getParam<MooseEnum>("method").getEnum<Method>()),
    _aitken_formula(getParam<MooseEnum>("aitken_formula").getEnum<AitkenFormula>()),
    _x_prev(_sys.addVector(name() + "_secant_x_prev", false, libMesh::PARALLEL)),
    _f_prev(_sys.addVector(name() + "_secant_f_prev", false, libMesh::PARALLEL)),
    _num_its(0),
    _aitken_factor(_initial_relaxation_factor)
{
  for (const auto & var_name : _var_names)
    if (!_fe_problem.hasSolverVariable(var_name))
      paramError("variables", "Variable '", var_name, "' is not a solver variable.");
    else if (&_fe_problem.getVariable(0, var_name).sys() != &_sys)
      paramError("variables", "All variables must belong to the same solver system.");

  if (_method != Method::AITKEN && isParamSetByUser("aitken_formula"))
    paramWarning("aitken_formula", "This parameter is only used with 'method = aitken'.");
  if (_method == Method::ANDERSON && isParamSetByUser("initial_relaxation_factor"))
    paramWarning("initial_relaxation_factor",
                 "This parameter is not used with 'method = anderson', use 'anderson_mixing'.");
  if (_method != Method::ANDERSON && isParamSetByUser("anderson_mixing"))
    paramWarning("anderson_mixing", "This parameter is only used with 'method = anderson'.");

  // The iterate entering the current multi-system fixed point iteration is needed
  _fe_problem.needsPreviousMultiSystemFixedPointIterationSolution(true, _sys.number());
}

void
MultiSystemSecantTransformation::initialSetup()
{
  gatherDofs();
}

void
MultiSystemSecantTransformation::timestepSetup()
{
  // A failed solve does not reach MULTISYSTEM_FIXED_POINT_END, so the history must also be
  // discarded before a (possibly repeated) time step
  _num_its = 0;
}

void
MultiSystemSecantTransformation::meshChanged()
{
  _num_its = 0;
  gatherDofs();
}

void
MultiSystemSecantTransformation::gatherDofs()
{
  std::vector<std::string> var_names(_var_names.begin(), _var_names.end());
  AllLocalDofIndicesThread aldit(_fe_problem, var_names);
  const libMesh::ConstElemRange & elem_range = *_fe_problem.mesh().getActiveLocalElementRange();
  Threads::parallel_reduce(elem_range, aldit, _fe_problem.numThreads());
  const auto & dofs = aldit.getDofIndices();
  _dofs.assign(dofs.begin(), dofs.end());
}

void
MultiSystemSecantTransformation::execute()
{
  if (_fe_problem.getCurrentExecuteOnFlag() == EXEC_MULTISYSTEM_FIXED_POINT_END)
  {
    _num_its = 0;
    return;
  }

  // The flag is executed after the solve of every system: only act after the solve of ours
  if (_fe_problem.justSolvedSolverSys() != _sys.number())
    return;

  auto & solution = _sys.solution();
  // Iterate x_k that entered this iteration, and its evaluation f(x_k) from this iteration
  const auto & x = _sys.solutionState(1, Moose::SolutionIterationType::MultiSystemFixedPoint);
  const auto & f = solution;

  std::vector<Number> x_k(_dofs.size()), f_k(_dofs.size()), x_km1(_dofs.size()),
      f_km1(_dofs.size());
  x.get(_dofs, x_k);
  f.get(_dofs, f_k);
  _x_prev.get(_dofs, x_km1);
  _f_prev.get(_dofs, f_km1);

  std::vector<Number> new_values(_dofs.size());
  // Coefficients of the update, for output
  std::ostringstream coefficients;
  coefficients << std::setprecision(8);
  if (_method == Method::PER_DOF)
  {
    Real step_min = std::numeric_limits<Real>::max();
    Real step_max = std::numeric_limits<Real>::lowest();
    for (const auto i : index_range(_dofs))
    {
      // Scalar secant step x_{k+1} = x_k + s r_k with s = -dx / dr. Two evaluations are needed
      // for the slope, so the first iteration is a Picard step with w_0. A null dr falls back
      // to a plain Picard step (s = 1)
      const Real r_k = f_k[i] - x_k[i];
      Real step = _num_its == 0 ? _initial_relaxation_factor : 1.0;
      if (_num_its > 0)
      {
        const Real dr = r_k - (f_km1[i] - x_km1[i]);
        if (!MooseUtils::absoluteFuzzyEqual(dr, 0))
          step = -(x_k[i] - x_km1[i]) / dr;
      }
      new_values[i] = x_k[i] + step * r_k;
      step_min = std::min(step_min, step);
      step_max = std::max(step_max, step);
    }
    _communicator.min(step_min);
    _communicator.max(step_max);
    coefficients << "per-DOF relaxation factor min = " << step_min << ", max = " << step_max;
  }
  else
  {
    // Global (all DOFs) scalar products needed by the vector secant methods, with
    // r = f(x) - x, dx = x_k - x_{k-1}, dr = r_k - r_{k-1}
    Real dr_dr = 0, dx_dr = 0, rkm1_dr = 0, rk_dr = 0, scale = 0;
    if (_num_its > 0)
    {
      for (const auto i : index_range(_dofs))
      {
        const Real r_k = f_k[i] - x_k[i];
        const Real r_km1 = f_km1[i] - x_km1[i];
        const Real dr = r_k - r_km1;
        const Real dx = x_k[i] - x_km1[i];
        dr_dr += dr * dr;
        dx_dr += dx * dr;
        rkm1_dr += r_km1 * dr;
        rk_dr += r_k * dr;
        scale += r_k * r_k + r_km1 * r_km1;
      }
      std::vector<Real> sums = {dr_dr, dx_dr, rkm1_dr, rk_dr, scale};
      _communicator.sum(sums);
      dr_dr = sums[0];
      dx_dr = sums[1];
      rkm1_dr = sums[2];
      rk_dr = sums[3];
      scale = sums[4];
    }
    // The residual difference is considered null when it is at round-off level
    const Real eps = std::numeric_limits<Real>::epsilon();
    const bool dr_ok = _num_its > 0 && dr_dr > eps * eps * scale;

    if (_method == Method::AITKEN)
    {
      if (_num_its == 0)
        _aitken_factor = _initial_relaxation_factor;
      else if (dr_ok)
      {
        // The previous factor is kept if dr is at round-off level
        if (_aitken_formula == AitkenFormula::RECURSIVE)
          _aitken_factor = -_aitken_factor * rkm1_dr / dr_dr;
        else
          _aitken_factor = -dx_dr / dr_dr;
      }

      for (const auto i : index_range(_dofs))
        new_values[i] = x_k[i] + _aitken_factor * (f_k[i] - x_k[i]);
      coefficients << "Aitken w_k = " << _aitken_factor;
    }
    else
    {
      // Anderson acceleration with m = 1 and mixing parameter beta:
      // gamma = argmin |r_k - gamma dr| = r_k.dr / |dr|^2
      // x_{k+1} = x_k + beta r_k - gamma (dx + beta dr)
      // The first iteration (AA(0)) and a degenerate dr (gamma = 0) are Picard steps with beta
      const Real beta = _anderson_mixing;
      const Real gamma = dr_ok ? rk_dr / dr_dr : 0.0;
      for (const auto i : index_range(_dofs))
      {
        const Real r_k = f_k[i] - x_k[i];
        const Real dr = r_k - (f_km1[i] - x_km1[i]);
        const Real dx = x_k[i] - x_km1[i];
        new_values[i] = x_k[i] + beta * r_k - gamma * (dx + beta * dr);
      }
      coefficients << "Anderson beta = " << beta << ", gamma = " << gamma;
    }
  }

  if (_verbose)
  {
    // Effective relaxation factor: projection of the update on the Picard direction r_k. It
    // equals the relaxation factor w_k whenever the update is parallel to r_k, e.g. a scalar or
    // pseudo-scalar problem where the Jacobian of f is a multiple of the identity
    Real update_dot_r = 0, r_dot_r = 0;
    for (const auto i : index_range(_dofs))
    {
      const Real r_k = f_k[i] - x_k[i];
      update_dot_r += (new_values[i] - x_k[i]) * r_k;
      r_dot_r += r_k * r_k;
    }
    std::vector<Real> sums = {update_dot_r, r_dot_r};
    _communicator.sum(sums);
    _console << name() << ": iteration " << _num_its << ", " << coefficients.str()
             << ", w_eff = "
             << (sums[1] > 0 ? std::to_string(sums[0] / sums[1]) : std::string("n/a")) << std::endl;
  }

  _x_prev.insert(x_k, _dofs);
  _f_prev.insert(f_k, _dofs);
  _x_prev.close();
  _f_prev.close();

  solution.insert(new_values, _dofs);
  solution.close();
  _sys.update();

  if (auto * const linear_sys = dynamic_cast<LinearSystem *>(&_sys))
    linear_sys->computeGradients();

  ++_num_its;
}
