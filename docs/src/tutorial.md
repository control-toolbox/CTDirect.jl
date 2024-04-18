# [First example: Goddard problem](@id goddard)

We start with the well-known so-called [Goddard](http://en.wikipedia.org/wiki/Robert_H._Goddard) problem, which models the ascent of a rocket through the atmosphere (see for instance [^1],[^2]).

[^1]: R.H. Goddard. A Method of Reaching Extreme Altitudes, volume 71(2) of Smithsonian Miscellaneous Collections. Smithsonian institution, City of Washington, 1919.

[^2]: H. Seywald and E.M. Cliff. Goddard problem in presence of a dynamic pressure limit. Journal of Guidance, Control, and Dynamics, 16(4):776–781, 1993.

```@raw html
<img src="./assets/goddard.png" style="display: block; margin: 0 auto 20px auto>
```

We restrict here ourselves to  vertical (monodimensional) trajectories, and the state variables are the altitude, speed and mass of the rocket during the flight, for a total dimension of 3. The rocket is subject to gravity, thrust and drag forces. The final time is free, and the objective is here to reach a maximal altitude with a given fuel consumption, i.e a fixed final mass. We also impose a speed limit. All units are renormalized.

## Basic solve

First import the CTDirect and CTBase modules
```@example main
using CTDirect
using CTBase
```

Then define the OCP for the Goddard problem. Note that the free final time is modeled as an *optimization variable*, and has both a lower bound to prevent the optimization getting stuck at tf=0. In this particular case an upper bound is not needed for the final time since the final mass is prescribed.
**+++ can we get here the nice table with the green checkmarks for the ocp ? it seems to appear in the terminal instead**
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
@def ocp1 begin
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

We can solve the problem directly using the default options, i.e 100 time steps and a constant initial guess set at 0.1. 
```@example main
sol1 = solve(ocp1, print_level=5)
nothing # hide
```
Then plot the solution with the state and control variables, as well as the costate recovered from the lagrange multipliers of the discretized problem. 
```@example main
plot(sol1)
```

## Initial guess options

+++ variantes initial guess: constant
sol2
, fonctionnel
sol3
, mixed
sol4

We can also use a so-called *warmstart* strategy and use an existing solution as initial guess (note that the OCP solution returned by the **solve** call is functional, thus it is not necessary to use the same time grid).
```@example main
sol4 = solve(ocp1, grid_size=200, print_level=5, init=OptimalControlInit(sol1))
nothing # hide
```
```@example main
plot(sol4)
```

## Working explicitely on the discretized problem
+++ variantes appel: transcription puis solve docp, avec init changee dans le docp (ex: mixed) et passee au solve (ex: warmstart)

