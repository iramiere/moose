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
 * Applies a secant update to selected variables at the end of each multi-system fixed point
 * iteration, either independently for each degree of freedom or with Aitken's dynamic relaxation.
 */
class MultiSystemSecantTransformation : public GeneralUserObject
{
public:
  static InputParameters validParams();

  MultiSystemSecantTransformation(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void timestepSetup() override;
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
  /// Relaxation factor applied to the per-DOF secant update, or initial Aitken relaxation factor
  const Real _relaxation_factor;

  /// Secant update method
  enum class Method
  {
    PER_DOF,
    AITKEN
  };
  const Method _method;

  /// Iterate x_{k-1} entering the previous multi-system fixed point iteration
  NumericVector<Number> & _x_prev;
  /// Evaluation f(x_{k-1}) of the previous multi-system fixed point iteration
  NumericVector<Number> & _f_prev;

  /// Number of multi-system fixed point iterations recorded since the last reset
  unsigned int _num_its;
  /// Current Aitken relaxation factor
  Real _aitken_factor;

  /// Local degrees of freedom of the transformed variables
  std::vector<dof_id_type> _dofs;
};
