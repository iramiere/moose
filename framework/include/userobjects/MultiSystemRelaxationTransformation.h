//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"

/**
 * Applies a static relaxation to selected variables right after the solve of their system in each
 * multi-system fixed point iteration (execute_on MULTISYSTEM_FIXED_POINT_TRANSFORMATION, executed
 * by the solve loop after every system solve):
 * x_{k+1} = x_k + w (f(x_k) - x_k), with x_k the iterate entering iteration k and f(x_k) its
 * evaluation.
 */
class MultiSystemRelaxationTransformation : public GeneralUserObject
{
public:
  static InputParameters validParams();

  MultiSystemRelaxationTransformation(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void meshChanged() override;

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}

protected:
  /// Gather the local degrees of freedom of the transformed variables
  void gatherDofs();

  /// Solver system holding the transformed variables
  SystemBase & _sys;
  /// Names of the transformed variables
  const std::vector<VariableName> & _var_names;
  /// Relaxation factor w
  const Real _relaxation_factor;

  /// Local degrees of freedom of the transformed variables
  std::vector<dof_id_type> _dofs;
};
