# [First example: Goddard problem](@id goddard)

We start with the well-known so-called [Goddard](http//en.wikipedia.org/wiki/Robert_H._Goddard) problem, which models the ascent of a rocket through the atmosphere (see for instance [1],[2]).

+++pic goddard.jpg

We restrict here ourselves to  vertical (monodimensional) trajectories, and the state variables are the altitude, speed and mass of the rocket during the flight, for a total dimension of 3. The rocket is subject to gravity, thrust and drag forces. The final time is free, and the objective is here to reach a maximal altitude with a given fuel consumption. All units are renormalized.

+++ OCP latex
```math

```

First import the module
```@setup main
using Plots
using CTDirect
```
Then define the OCP for the Goddard problem. Note that the free final time is modeled as an optimization variable, hence the argument (variable=true).
```@example main
ocp1 = Model(variable=true)
Cd = 310
Tmax = 3.5
β = 500
b = 2
r0 = 1
v0 = 0
vmax = 0.1
m0 = 1
mf = 0.6
x0 = [ r0, v0, m0 ]
state!(ocp1, 3)
control!(ocp1, 1)
variable!(ocp1, 1)
time!(ocp1, 0, Index(1))
constraint!(ocp1, :initial, x0, :initial_constraint)
constraint!(ocp1, :final, Index(3), mf, :final_constraint)
constraint!(ocp1, :control, u->u, 0, 1, :control_bounds)
constraint!(ocp1, :state, 1:2, [r0,v0,mf], [Inf, vmax, m0], :state_bounds)
constraint!(ocp1, :variable, Index(1), 0.01, 10, :variable_tf_bounds)
objective!(ocp1, :mayer,  (x0, xf, v) -> xf[1], :max)
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
dynamics!(ocp, (x, u, v) -> F0(x) + u*F1(x) )
```

We can then solve the problem directly
```@example main
sol0 = solve(ocp, print_level=0)
```
and plot the solution
```@example main
plot(sol0)
```
+++ abstract ocp

+++ variantes initial guess: constant, fonctionnel, mixed, warmstart

+++ variantes appel: transcription puis solve docp, avec init changee dans le docp (ex: mixed) et passee au solve (ex: warmstart)



References

[1] R.H. Goddard. A Method of Reaching Extreme Altitudes, volume 71(2) of Smithsonian Miscellaneous Collections. Smithsonian institution, City of Washington, 1919.

[2] H. Seywald and E.M. Cliff. Goddard problem in presence of a dynamic pressure limit. Journal of Guidance, Control, and Dynamics, 16(4):776–781, 1993.