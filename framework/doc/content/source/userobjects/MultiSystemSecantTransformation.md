# MultiSystemSecantTransformation

!syntax description /UserObjects/MultiSystemSecantTransformation

## Description

This user object accelerates the multi-system fixed point iteration of [FEProblemSolve.md] by
transforming the [!param](/UserObjects/MultiSystemSecantTransformation/variables) at the end of each
iteration. Let $x_k$ be the values entering iteration $k$, $f(x_k)$ the values obtained after solving every
system in iteration $k$, and $r_k = f(x_k) - x_k$ the fixed point residual. Two secant updates are
available through [!param](/UserObjects/MultiSystemSecantTransformation/method).

With `per_dof`, each degree of freedom is updated independently, as in [SecantSolve.md]:

!equation
x_{k+1} = x_k - \dfrac{r_k (x_k - x_{k-1})}{r_k - r_{k-1}}

which is then relaxed with the [!param](/UserObjects/MultiSystemSecantTransformation/relaxation_factor)
$\omega$ as $\omega x_{k+1} + (1 - \omega) x_k$. Degrees of freedom whose secant denominator vanishes keep
their Picard update.

With `aitken`, a single relaxation factor shared by all transformed degrees of freedom is computed with
Aitken's delta-squared method:

!equation
\omega_k = -\omega_{k-1} \dfrac{r_{k-1} \cdot (r_k - r_{k-1})}{\lVert r_k - r_{k-1} \rVert^2}, \qquad
x_{k+1} = x_k + \omega_k r_k

where the initial factor $\omega_0$ is the
[!param](/UserObjects/MultiSystemSecantTransformation/relaxation_factor). Because the factor is computed
from inner products over all the transformed degrees of freedom, this method is better suited than
`per_dof` to problems where the degrees of freedom are strongly coupled, for example by diffusion.

With `per_dof`, the first iteration of each multi-system fixed point loop is a plain Picard update,
because the secant slope requires two evaluations.

The transformation is executed on `MULTISYSTEM_FIXED_POINT_ITERATION_END`, and the iteration history is
discarded on `MULTISYSTEM_FIXED_POINT_END`, at each time step, and when the mesh changes. All variables
must belong to the same solver system. The secant update is computed from the solution after all the
systems have been solved in an iteration, so systems solved after the one holding the transformed
variables see the untransformed values within that iteration.

## Example Input Syntax

!listing test/tests/multisystem/picard/nonlinearfe_nonlinearfe/multi_system_secant.i block=UserObjects

!syntax parameters /UserObjects/MultiSystemSecantTransformation

!syntax inputs /UserObjects/MultiSystemSecantTransformation

!syntax children /UserObjects/MultiSystemSecantTransformation
