# This example aims to replicate the test cas "Single-block thermomechanics example" with two ways coupling
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
  displacements = 'disp_x disp_y disp_z'
  # type = FileMesh
  # file = single_block_mesh_refined.e
[]

[Variables]
  # We solve for the temperature and the displacements
  [T]
    initial_condition = 100
  []
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]


[Kernels]
  [htcond] #Heat conduction equation
    type = HeatConduction
    variable = T
    use_displaced_mesh = true
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
    type= ComputeThermalExpansionEigenstrain
    #thermal_expansion_coeff = 1e-5
    thermal_expansion_coeff = 1e-3
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
  #  automatic_scaling =true
  #  scaling_group_variables ='disp_x disp_y disp_z'
  #  off_diagonals_in_auto_scaling= true
  # resid_vs_jac_scaling_param=1
  # line_search= none
  verbose=true

   petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
   petsc_options_value = 'hypre boomeramg 101'
   ### for Newton without preconditioning
  # petsc_options_iname = '-pc_type  -ksp_gmres_restart'
  # petsc_options_value = 'none  101'
   ### for tuning e_rel
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_err'
  # petsc_options_value = 'hypre boomeramg 101 1'

  # for conditioning analysis (with Newton solve)
  # petsc_options = '-pc_svd_monitor'
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'svd'

  l_max_its = 30
  nl_max_its = 10
  nl_rel_tol = 1e-08
  l_abs_tol = 1e-15
  l_tol = 1e-12
[]

[Outputs]
  exodus = true
  perf_graph = true
[]
