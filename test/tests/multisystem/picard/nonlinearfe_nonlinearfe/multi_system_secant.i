# Solves problem in a multi-system fixed point solve, accelerating the iteration with
# Aitken's dynamic relaxation of v at the end of each multi-system fixed point iteration.

!include multi_system.i

[UserObjects]
  [secant]
    type = MultiSystemSecantTransformation
    variables = 'v'
    method = aitken
  []
[]
