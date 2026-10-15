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

registerMooseObject("MooseApp", MultiSystemSecantTransformation);

InputParameters
MultiSystemSecantTransformation::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Applies a secant update to selected variables at the end of each multi-system fixed point "
      "iteration, either independently for each degree of freedom or with Aitken's dynamic "
      "relaxation.");
  params.addRequiredParam<std::vector<VariableName>>(
      "variables", "Solver variables to transform. They must all belong to the same solver system.");
  params.addRangeCheckedParam<Real>(
      "relaxation_factor",
      1.0,
      "relaxation_factor>0 & relaxation_factor<2",
      "Relaxation factor applied to the secant update with 'per_dof', or initial relaxation "
      "factor of the iteration with 'aitken'");
  MooseEnum method("per_dof aitken", "per_dof");
  params.addParam<MooseEnum>(
      "method",
      method,
      "Secant update to use: 'per_dof' applies an independent secant update to each degree of "
      "freedom, 'aitken' computes a single dynamic relaxation factor for all degrees of freedom "
      "using Aitken's delta-squared method");

  // The transformation is applied at the end of each iteration, and the iteration history is
  // reset at the end of the loop
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum = {EXEC_MULTISYSTEM_FIXED_POINT_ITERATION_END, EXEC_MULTISYSTEM_FIXED_POINT_END};
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

MultiSystemSecantTransformation::MultiSystemSecantTransformation(
    const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _sys(_fe_problem.getVariable(0, getParam<std::vector<VariableName>>("variables")[0]).sys()),
    _var_names(getParam<std::vector<VariableName>>("variables")),
    _relaxation_factor(getParam<Real>("relaxation_factor")),
    _method(getParam<MooseEnum>("method").getEnum<Method>()),
    _x_prev(_sys.addVector(name() + "_secant_x_prev", false, libMesh::PARALLEL)),
    _f_prev(_sys.addVector(name() + "_secant_f_prev", false, libMesh::PARALLEL)),
    _num_its(0),
    _aitken_factor(_relaxation_factor)
{
  for (const auto & var_name : _var_names)
    if (!_fe_problem.hasSolverVariable(var_name))
      paramError("variables", "Variable '", var_name, "' is not a solver variable.");
    else if (&_fe_problem.getVariable(0, var_name).sys() != &_sys)
      paramError("variables", "All variables must belong to the same solver system.");

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
  if (_method == Method::PER_DOF)
    for (const auto i : index_range(_dofs))
    {
      // Two evaluations are needed for the secant slope, so the first iteration is a Picard
      // update
      Real new_value = f_k[i];
      if (_num_its > 0)
      {
        // Avoid a zero denominator
        const Real denominator = f_k[i] - x_k[i] - f_km1[i] + x_km1[i];
        if (!MooseUtils::absoluteFuzzyEqual(denominator, 0))
          new_value = x_k[i] - (f_k[i] - x_k[i]) * (x_k[i] - x_km1[i]) / denominator;
        new_value = _relaxation_factor * new_value + (1 - _relaxation_factor) * x_k[i];
      }
      new_values[i] = new_value;
    }
  else
  {
    // Aitken's dynamic relaxation factor:
    // w_k = -w_{k-1} r_{k-1}.(r_k - r_{k-1}) / |r_k - r_{k-1}|^2 with residuals r = f(x) - x
    if (_num_its == 0)
      _aitken_factor = _relaxation_factor;
    else
    {
      Real numerator = 0;
      Real denominator = 0;
      for (const auto i : index_range(_dofs))
      {
        const Real r_km1 = f_km1[i] - x_km1[i];
        const Real delta_r = f_k[i] - x_k[i] - r_km1;
        numerator += r_km1 * delta_r;
        denominator += delta_r * delta_r;
      }
      _communicator.sum(numerator);
      _communicator.sum(denominator);
      // Keep the previous factor if the residual did not change
      if (!MooseUtils::absoluteFuzzyEqual(denominator, 0))
        _aitken_factor = -_aitken_factor * numerator / denominator;
    }

    for (const auto i : index_range(_dofs))
      new_values[i] = x_k[i] + _aitken_factor * (f_k[i] - x_k[i]);
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
