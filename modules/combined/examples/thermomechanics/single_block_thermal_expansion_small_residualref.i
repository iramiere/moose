# This example aims to replicate the test cas "Single-block thermomechanics example"
# of the paper 'Evaluation of coupling approaches for thermomechanical simulations'
# S.R. Novascone∗, B.W. Spencer, J.D. Hales, R.L. Williamson
# NED 2015

[Problem]
  type = ReferenceResidualProblem
  reference_vector = 'ref'
  extra_tag_vectors = 'ref'
[]


[Mesh]
  #rectangular block
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
    extra_vector_tags = 'ref'
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
    extra_vector_tags = 'ref'
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

[Executioner]
  type = Transient
  num_steps = 1
  solve_type = PJFNK
  #automatic_scaling =true
  #scaling_group_variables ='disp_x disp_y disp_z'
  #off_diagonals_in_auto_scaling= true
  #resid_vs_jac_scaling_param=1
  verbose=true

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 101'
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_err'
  #petsc_options_value = 'hypre boomeramg 101 1'
  l_max_its = 30
  nl_max_its = 10
  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10
  l_tol = 1e-12
[]

[Outputs]
  exodus = true
  perf_graph = true
[]
