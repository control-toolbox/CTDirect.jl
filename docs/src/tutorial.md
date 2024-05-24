# First example: Goddard problem

![R.H. Goddard](./assets/goddard.png)

We start with the well-known so-called [Goddard](http://en.wikipedia.org/wiki/Robert_H._Goddard) problem, which models the ascent of a rocket through the atmosphere (see for instance [^1],[^2]).
We restrict here ourselves to  vertical (monodimensional) trajectories, and the state variables are the altitude, speed and mass of the rocket during the flight, for a total dimension of 3. The rocket is subject to gravity, thrust and drag forces. The final time is free, and the objective is here to reach a maximal altitude with a given fuel consumption, i.e a fixed final mass. We also impose a speed limit. All units are renormalized.

[^1]: R.H. Goddard. A Method of Reaching Extreme Altitudes, volume 71(2) of Smithsonian Miscellaneous Collections. Smithsonian institution, City of Washington, 1919.

[^2]: H. Seywald and E.M. Cliff. Goddard problem in presence of a dynamic pressure limit. Journal of Guidance, Control, and Dynamics, 16(4):776–781, 1993.

## Problem definition

First import the CTDirect and CTBase modules
```@example main
using CTDirect
using CTBase
```

Then define the OCP for the Goddard problem. Note that the free final time is modeled as an *optimization variable*, and has both a lower bound to prevent the optimization getting stuck at tf=0. In this particular case an upper bound is not needed for the final time since the final mass is prescribed.
```@example main
Cd = 310
β = 500
Tmax = 3.5
b = 2
vmax = 0.1
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
@def ocp begin
    tf ∈ R, variable
    t ∈ [ 0, tf ], time
    x ∈ R^3, state
    u ∈ R, control
    0.1 ≤ tf ≤ Inf
    r = x[1]
    v = x[2]
    m = x[3]
    x(0) == [1, 0, 1]
    m(tf) == 0.6
    1 ≤ r(t) ≤ 1.1
    0 ≤ v(t) ≤ vmax
    0 ≤ u(t) ≤ 1
    ẋ(t) == F0(x(t)) + u(t)*F1(x(t))
    r(tf) → max
end
nothing # hide
```

## Basic solve

We can solve the problem directly using the default options.
```@example main
sol1 = solve(ocp, print_level=5)
nothing # hide
```
Then plot the solution with the state and control variables, as well as the costate recovered from the lagrange multipliers of the discretized problem. 
```@example main
plot(sol1)
```

The most common option for **solve** is the number of time steps for the discretized problem (default 100), that can be set with the argument *grid_size*. A larger grid size will increase the computational cost, while a smaller value may lead to a very coarse solution.

## Initial guess options

The function **solve** uses a default constant initialisation of 0.1 for all variables. More advanced options include constant and/or functional initialisation for each individual state or control component, as well as reusing an existing solution, also known as *warm start*[^3].

[^3]: Currently only the primal variables are reused for the warm start, not the lagrange multipliers. It should be noted that previous experiments with the Bocop software seemed to indicate that initializing also the multipliers gave little benefit.

Let us start with the simplest case, constant initialisation.
```@example main
x_const = [1.05, 0.2, 0.8]
u_const = 0.5
v_const = 0.15
init1 = OCPInit(state=x_const, control=u_const, variable=v_const)
sol2 = solve(ocp, print_level=0, init=init1)
println("Objective ", sol2.objective, " after ", sol2.iterations, " iterations")
``` 

Now we illustrate the functional initialisation, with some random functions. Note that we only consider the state and control variables, since the optimization variables are scalar and therefore a functional initialisation is not relevant. In the example notice that the call to **OCPInit** does not provide an argument for the optimization variables, therefore the default initial guess will be used.  
```@example main
x_func = t->[1+t^2, sqrt(t), 1-t]
u_func = t->(cos(t)+1)*0.5
init2 = OCPInit(state=x_func, control=u_func)
sol3 = solve(ocp, print_level=0, init=init2)
println("Objective ", sol3.objective, " after ", sol3.iterations, " iterations")
``` 
More generally, the default, constant and functional initialisations can be mixed, as shown in the example below that uses a functional initial guess for the state, a constant initial guess for the control, and the default initial guess for the optimization variables. 
```@example main
init3 = OCPInit(state=x_func, control=u_const)
sol4 = solve(ocp, print_level=0, init=init3)
println("Objective ", sol4.objective, " after ", sol4.iterations, " iterations")
``` 

Finally, we can also use a so-called *warmstart* strategy and use an existing solution as initial guess (note that the OCP solution returned by the **solve** call is functional, thus it is not necessary to use the same time grid). Notice that the objective and constraint violation values start much closer to the solution than with the previous initialisations.
```@example main
sol4 = solve(ocp, grid_size=200, print_level=5, init=sol1)
nothing # hide
```
```@example main
plot(sol4)
```

## The discretized problem

Instead of calling **solve** directly on the OCP problem, you can first obtain the discretized problem (DOCP) by calling **directTranscription**, then call **solve** on the DOCP. The resulting solution of the discretized problem can be used to generate the corresponding OCP solution with **OCPSolutionFromDOCP**.
```@example main
docp = directTranscription(ocp, grid_size=100)
dsol = solve(docp, print_level=5)
sol5 = OCPSolutionFromDOCP(docp, dsol)
nothing # hide
```
The initial guess can be passed to **solve** same as before.
```@example main
dsol = solve(docp, print_level=0, init=sol1)
sol6 = OCPSolutionFromDOCP(docp, dsol)
println("Objective ", sol6.objective, " after ", sol6.iterations, " iterations")
```
Another possibility is to set the initial guess associated to the DOCP, using the function **setDOCPInit**.
```@example main
setDOCPInit(docp, sol1)
dsol = solve(docp, print_level=5)
nothing # hide
```
Finally, the direct transcription also accept an initial guess.
```@example main
docp = directTranscription(ocp, grid_size=100, init=sol1)
dsol = solve(docp, print_level=5)
nothing # hide
```
