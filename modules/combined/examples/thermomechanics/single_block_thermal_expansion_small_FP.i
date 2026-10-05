# This example aims to replicate the test cas "Single-block thermomechanics example" with one way coupling
# of the paper 'Evaluation of coupling approaches for thermomechanical simulations'
# S.R. Novascone∗, B.W. Spencer, J.D. Hales, R.L. Williamson
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
  # displacements = 'disp_x disp_y disp_z'
  # type = FileMesh
  # file = single_block_mesh.e
[]

[Problem]
  nl_sys_names = 'temp disp'
[]

[Variables]
  # We solve for the temperature and the displacements
  [T]
    initial_condition = 100
    scaling = '1'
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
    # type = HeatConduction
    type = MatDiffusion
    diffusivity = 'thermal_conductivity'
    variable = T
    use_displaced_mesh = false
  []
  # [./TensorMechanics] #Action that creates equations for disp_x and disp_y
  #   displacements = 'disp_x disp_y'
  # [../]
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
  [all_T] #Temperature on the boundary is fixed to 1100 K
    type = FunctionDirichletBC
    variable = T
    boundary = 'front back top bottom left right'
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
    type= ComputeThermalExpansionEigenstrain
    thermal_expansion_coeff = 1e-5
    temperature = T
    stress_free_temperature = 100
    eigenstrain_name = thermal_expansion
  []
[]

# [Preconditioning]
#   [smp]
#     type = SMP
#     full = true
#   []
# []

[Executioner]
  type = Transient
  num_steps = 1
  abort_on_solve_fail = true
  solve_type = NEWTON

  line_search = 'none'
  verbose=true


  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = 'multi_sys'

  nonlinear_convergence = 'temp_conv disp_conv'

  # nl_max_its = 10
  # nl_rel_tol = 1e-08
  # nl_abs_tol = 1e-10

  l_max_its = 80
  l_abs_tol = 1e-15
  l_tol = 1e-12

[]

[Preconditioning]
  [temp_param]
    type = SMP
    # petsc_options = '-ksp_view'
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type'
    petsc_options_value = 'hypre   boomeramg  cg'
    # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type -mat_mffd_err'
    # petsc_options_value = 'hypre   boomeramg  cg 1'
    nl_sys = 'temp'
  []
  [disp_param]
    type = SMP
    # petsc_options = '-ksp_view'
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type'
    petsc_options_value = 'hypre   boomeramg  cg'
    # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type -mat_mffd_err'
    # petsc_options_value = 'hypre   boomeramg  cg 1'
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
    nl_abs_tol = 1e-10
    nl_rel_tol = 1e-8
    nl_max_its = 10
  []

  [multi_sys]
    # type = DefaultNonlinearConvergence
    # nl_max_its = 10
    # nl_rel_tol = 1e-8
    type = ParsedConvergence
    # convergence_expression = '(nl_temp + nl_disp) < 1e-8'
    convergence_expression = '(nl_temp < 1e-8) & (nl_disp < 1e-8)'
    symbol_names = 'nl_temp nl_disp'
    symbol_values = 'nl_temp nl_disp'
  []
[]

[Postprocessors]
  [nl_temp]
    type = Residual
    residual_type = 'compute'
    solver_sys = 'temp'
    execute_on = 'MULTISYSTEM_FIXED_POINT_ITERATION_END'
  []
  [nl_disp]
    type = Residual
    residual_type = 'compute'
    solver_sys  = 'disp'
    execute_on = 'MULTISYSTEM_FIXED_POINT_ITERATION_END'
  []
[]

[Outputs]
  exodus = true
  perf_graph = true
  # color = false
[]
