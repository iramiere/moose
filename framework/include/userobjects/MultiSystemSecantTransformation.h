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
 * Applies a secant-type acceleration to selected variables right after the solve of their system
 * in each multi-system fixed point iteration (execute_on MULTISYSTEM_FIXED_POINT_TRANSFORMATION,
 * executed by the solve loop after every system solve). With x_k the iterate of the system entering
 * iteration k, f(x_k) its evaluation and r_k = f(x_k) - x_k the fixed point residual, the available
 * updates are:
 *  - per_dof:   independent scalar secant (Steffensen-like) update of each degree of freedom
 *  - aitken:    x_{k+1} = x_k + w_k r_k with a single scalar w_k for all the degrees of freedom
 *  - anderson:  Anderson acceleration with depth m = 1 (vector secant, residual minimization)
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
  /// Relaxation factor w_0 of the first iteration of per_dof and aitken (no secant information yet)
  const Real _initial_relaxation_factor;
  /// Mixing parameter beta of the Anderson acceleration, applied at every iteration
  const Real _anderson_mixing;
  /// Whether to print the coefficients of the update at each iteration
  const bool _verbose;

  /// Acceleration method
  enum class Method
  {
    PER_DOF,
    AITKEN,
    ANDERSON
  };
  const Method _method;

  /// Formula used to compute the Aitken relaxation factor
  enum class AitkenFormula
  {
    RECURSIVE,
    CROSSED_SECANT
  };
  const AitkenFormula _aitken_formula;

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
