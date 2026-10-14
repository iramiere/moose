# This example aims to replicate the test cas "Single-block thermomechanics example" with two ways coupling
# of the paper 'Evaluation of coupling approaches for thermomechanical simulations'
# S.R. Novascone, B.W. Spencer, J.D. Hales, R.L. Williamson
# NED 2015

[Mesh]
  # rectangular block
  [block]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = 5
    ymin = 0
    ymax = 1
    zmin = 0
    zmax = 1
    nx = 20
    ny = 4
    nz = 4
    elem_type = HEX8
  []
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  nl_sys_names = 'temp disp'
[]

[Variables]
  # We solve for the temperature and the displacements
  [T]
    initial_condition = 100
    solver_sys = 'temp'
  []
  [disp_x]
    solver_sys = 'disp'
  []
  [disp_y]
    solver_sys = 'disp'
  []
  [disp_z]
    solver_sys = 'disp'
  []
[]

[Kernels]
  [htcond] #Heat conduction equation
    type = HeatConduction
    variable = T
    use_displaced_mesh = true
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  displacements = 'disp_x disp_y disp_z'
  [block]
    strain = SMALL
    displacements = 'disp_x disp_y disp_z'
    eigenstrain_names = 'thermal_expansion'
    temperature = T
  []
[]

[Functions]
  [temperature_ramp]
    type = PiecewiseLinear
    x = '0 1'
    y = '100 1100'
  []
[]

[BCs]
  [fix_temp_right] #Temperature on the right boundary is fixed to 100 K
    type = DirichletBC
    variable = T
    boundary = right
    value = 100
  []
  [ramp_temp_left] #Temperature on the left boundary is fixed to 1100 K
    type = FunctionDirichletBC
    variable = T
    boundary = left
    function = temperature_ramp
  []
  [fix_disp_x] #Displacements in the x-direction are fixed at the left
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_disp_y] #Displacements in the y-direction are fixed on all the boundary (axial displavcement only)
    type = DirichletBC
    variable = disp_y
    boundary = 'front back top bottom left right'
    value = 0
  []
  [fix_disp_z] #Displacements in the z-direction are fixed on all the boundary (axial displavcement only)
    type = DirichletBC
    variable = disp_z
    boundary = 'front back top bottom left right'
    value = 0
  []
[]

[Materials]
  [thcond] #Thermal conductivity is set to 1 W/mK
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = 1
  []
  [elasticity_tensor] #Sets isotropic elastic constants
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.3
  []
  [stress] #We use linear elasticity
    type = ComputeLinearElasticStress
  []
  [thermal_strain]
    type = ComputeThermalExpansionEigenstrain
    thermal_expansion_coeff = 1e-5
    temperature = T
    stress_free_temperature = 100
    eigenstrain_name = thermal_expansion
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
  abort_on_solve_fail = true
  solve_type = NEWTON

  line_search = 'none'
  verbose = true

  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = 'multi_sys'

  nonlinear_convergence = 'temp_conv disp_conv'

  # reuse_preconditioner = true
  # reuse_preconditioner_max_linear_its = 25


  l_max_its = 80
  l_abs_tol = 1e-15
  l_tol = 1e-12
[]

[Preconditioning]
  [temp_param]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type'
    petsc_options_value = 'hypre   boomeramg  cg'
    nl_sys = 'temp'
  []
  [disp_param]
    type = SMP
    full = true
    petsc_options = '-options_left'
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type
    -snes_lag_jacobian -snes_lag_jacobian_persists'
    petsc_options_value = 'hypre   boomeramg  cg -2 true'
    nl_sys = 'disp'
  []
[]

[Convergence]
  [temp_conv]
    type = DefaultNonlinearConvergence
    # l_max_its = 80
    nl_max_its = 10
    nl_abs_tol = 1e-10
    nl_rel_tol = 1e-8
  []
  [disp_conv]
    type = DefaultNonlinearConvergence
    nl_abs_tol = 1e-8
    nl_rel_tol = 1e-8
    nl_max_its = 10
  []

  # for multi_system
  [it_conv]
    type = IterationCountConvergence
    max_iterations = 100
    verbose = true
  []
  [ms_rel_conv]
    type = MultiPostprocessorConvergence
    postprocessors = 'nl_temp_relative nl_disp_relative'
    tolerances = '1e-6 1e-6'
    # tolerances = '1e-12 1e-12'
    verbose = true
  []
  [ms_abs_conv]
    type = MultiPostprocessorConvergence
    postprocessors = 'nl_temp nl_disp'
    tolerances = '1e-8 1e-8'
    verbose = true
  []
  [multi_sys]
    type = ParsedConvergence
    convergence_expression = 'rel_conv | abs_conv'
    divergence_expression = 'it_conv'
    symbol_names = 'rel_conv abs_conv it_conv'
    symbol_values = 'ms_rel_conv ms_abs_conv it_conv'
    verbose = true
  []
[]

[Postprocessors]

  # residual
  [nl_temp_relative]
    type = FixedPointRelativeResidual
    solver_sys = 'temp'
    execute_on = 'MULTISYSTEM_FIXED_POINT_ITERATION_END'
  []
  [nl_disp_relative]
    type = FixedPointRelativeResidual
    solver_sys = 'disp'
    execute_on = 'MULTISYSTEM_FIXED_POINT_ITERATION_END'
  []
  [nl_temp]
    type = Residual
    residual_type = 'compute'
    solver_sys = 'temp'
    execute_on = 'MULTISYSTEM_FIXED_POINT_ITERATION_END'
  []
  [nl_disp]
    type = Residual
    residual_type = 'compute'
    solver_sys = 'disp'
    execute_on = 'MULTISYSTEM_FIXED_POINT_ITERATION_END'
  []

  [t_begin]
    type = PerfGraphData
    section_name = 'Root'
    data_type = TOTAL
    execute_on = 'TIMESTEP_BEGIN'
  []
  [t_total]
    type = PerfGraphData
    section_name = 'Root'
    data_type = TOTAL
    execute_on = 'MULTISYSTEM_FIXED_POINT_ITERATION_END'
  []
  [t_solve]
    type = PerfGraphData
    section_name = 'FEProblem::solve'
    data_type = TOTAL
    must_exist = false
    execute_on = 'MULTISYSTEM_FIXED_POINT_ITERATION_END'
  []
[]

[Outputs]
  exodus = true
  perf_graph = true
  perf_graph_live = true
  perf_graph_live_time_limit = 1.0
  print_linear_residuals = true
  console = true
  # [pgraph]
  #   type = PerfGraphOutput
  #   execute_on = 'final'
  #   level = 2
  #   heaviest_sections = 10
  # []
[]
