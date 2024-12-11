API description for Discretization methods

- from docp.jl
objective
- mayer cost: initial and final state
- lagrange state at final time

constraints
- set work array (can be empty if not needed)
- dynamics constraints (state and stage equations)
- path constraints (u, x, xu)
- point constraints (boundary and variable)

bounds
- NLP variables bounds: x, u, v
- NLP constraints bounds: dyn, path, point

initial guess
- set state/control at time step (assume constant regardless of actual parametrization)
- set variable (by convention at the end of NLP variables, so not disc dependent)

- from solution.jl
retrieve state, control, variable, multipliers
stage grid ?

- internal
get state at time step (scal/vect)
get control at time step (scal/vect)
get variable (scal/vect)
