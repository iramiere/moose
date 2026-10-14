//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FixedPointRelativeResidual.h"

#include "FEProblem.h"
#include "SubProblem.h"
#include "NonlinearSystem.h"

registerMooseObject("MooseApp", FixedPointRelativeResidual);

InputParameters
FixedPointRelativeResidual::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Reports the residual of 'solver_sys' recomputed at the current state of the solution "
      "(not the final residual of its last internal Newton solve), normalized by the initial "
      "residual of that system captured once at the first execution and held fixed "
      "thereafter. Intended as a relative convergence measure for multi-system fixed-point "
      "(Picard) iterations that properly reflects the state of the coupled solution, unlike "
      "nl_rel_tol which only measures whether the last internal Newton solve converged.");
  params.addParam<SolverSystemName>("solver_sys", "Name of the nonlinear system");
  return params;
}

FixedPointRelativeResidual::FixedPointRelativeResidual(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _reference_residual(0.0),
    _has_reference(false),
    _current_relative(0.0)
{
}

void
FixedPointRelativeResidual::execute()
{
  FEProblemBase * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);
  if (!fe_problem)
    mooseError("Dynamic cast to FEProblemBase failed in FixedPointRelativeResidual");

  auto & nl_sys = fe_problem->getNonlinearSystemBase(_sys.number());

  const Real current_residual = fe_problem->computeResidualL2Norm(nl_sys);

  if (!_has_reference)
  {
    _reference_residual = nl_sys.initialResidual();
    _has_reference = true;
  }

  _current_relative =
      (_reference_residual != 0.0) ? current_residual / _reference_residual : current_residual;
}

Real
FixedPointRelativeResidual::getValue() const
{
  return _current_relative;
}