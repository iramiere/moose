# Counts executions on the multi-system fixed point execute flags.

!include multi_system.i

[Postprocessors]
  [num_iteration_end]
    type = GeneralSetupInterfaceCount
    count_type = FINALIZE
    execute_on = 'MULTISYSTEM_FIXED_POINT_ITERATION_END'
  []
  [num_end]
    type = GeneralSetupInterfaceCount
    count_type = FINALIZE
    execute_on = 'MULTISYSTEM_FIXED_POINT_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    file_base = multi_system_end_csv
    show = 'num_iteration_end num_end'
    execute_on = 'FINAL'
  []
[]
