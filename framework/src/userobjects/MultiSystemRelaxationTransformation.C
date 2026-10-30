//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultiSystemRelaxationTransformation.h"
#include "AllLocalDofIndicesThread.h"
#include "FEProblemBase.h"
#include "LinearSystem.h"

#include "libmesh/numeric_vector.h"

registerMooseObject("MooseApp", MultiSystemRelaxationTransformation);

InputParameters
MultiSystemRelaxationTransformation::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Applies a static relaxation to selected variables right after the solve of their system "
      "in each multi-system fixed point iteration.");
  params.addRequiredParam<std::vector<VariableName>>(
      "variables", "Solver variables to relax. They must all belong to the same solver system.");
  params.addRequiredRangeCheckedParam<Real>(
      "relaxation_factor",
      "relaxation_factor>0 & relaxation_factor<2",
      "Relaxation factor w applied to the update, x + w (f(x) - x); 1 means no relaxation.");

  // The relaxation is applied right after the solve of the system holding the variables
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum = {EXEC_MULTISYSTEM_FIXED_POINT_TRANSFORMATION};
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

MultiSystemRelaxationTransformation::MultiSystemRelaxationTransformation(
    const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _sys(_fe_problem.getVariable(0, getParam<std::vector<VariableName>>("variables")[0]).sys()),
    _var_names(getParam<std::vector<VariableName>>("variables")),
    _relaxation_factor(getParam<Real>("relaxation_factor"))
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
MultiSystemRelaxationTransformation::initialSetup()
{
  gatherDofs();
}

void
MultiSystemRelaxationTransformation::meshChanged()
{
  gatherDofs();
}

void
MultiSystemRelaxationTransformation::gatherDofs()
{
  std::vector<std::string> var_names(_var_names.begin(), _var_names.end());
  AllLocalDofIndicesThread aldit(_fe_problem, var_names);
  const libMesh::ConstElemRange & elem_range = *_fe_problem.mesh().getActiveLocalElementRange();
  Threads::parallel_reduce(elem_range, aldit, _fe_problem.numThreads());
  const auto & dofs = aldit.getDofIndices();
  _dofs.assign(dofs.begin(), dofs.end());
}

void
MultiSystemRelaxationTransformation::execute()
{
  // The flag is executed after the solve of every system: only act after the solve of ours
  if (_fe_problem.justSolvedSolverSys() != _sys.number())
    return;

  auto & solution = _sys.solution();
  // Iterate x_k that entered this iteration, and its evaluation f(x_k) from this iteration
  const auto & x = _sys.solutionState(1, Moose::SolutionIterationType::MultiSystemFixedPoint);
  const auto & f = solution;

  std::vector<Number> x_k(_dofs.size()), f_k(_dofs.size()), new_values(_dofs.size());
  x.get(_dofs, x_k);
  f.get(_dofs, f_k);

  for (const auto i : index_range(_dofs))
    new_values[i] = x_k[i] + _relaxation_factor * (f_k[i] - x_k[i]);

  solution.insert(new_values, _dofs);
  solution.close();
  _sys.update();

  if (auto * const linear_sys = dynamic_cast<LinearSystem *>(&_sys))
    linear_sys->computeGradients();
}
