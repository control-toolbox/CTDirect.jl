# [First example: Goddard problem](@id goddard)

We start with the well-known so-called [Goddard](http://en.wikipedia.org/wiki/Robert_H._Goddard) problem, which models the ascent of a rocket through the atmosphere (see for instance [1],[2]).

[1] R.H. Goddard. A Method of Reaching Extreme Altitudes, volume 71(2) of Smithsonian Miscellaneous Collections. Smithsonian institution, City of Washington, 1919.

[2] H. Seywald and E.M. Cliff. Goddard problem in presence of a dynamic pressure limit. Journal of Guidance, Control, and Dynamics, 16(4):776–781, 1993.

```@raw html
<img src="./assets/goddard.png" style="display: block; margin: 0 auto 20px auto>
```

We restrict here ourselves to  vertical (monodimensional) trajectories, and the state variables are the altitude, speed and mass of the rocket during the flight, for a total dimension of 3. The rocket is subject to gravity, thrust and drag forces. The final time is free, and the objective is here to reach a maximal altitude with a given fuel consumption. All units are renormalized.

+++ OCP latex
```math

```

```@setup main
using Plots
using Plots.PlotMeasures
plot(args...; kwargs...) = Plots.plot(args...; kwargs..., leftmargin=25px)
```

First import the CTDirect module
```@example main
using CTDirect
```
Then define the OCP for the Goddard problem. Note that the free final time is modeled as an *optimization variable*.

```@example main
@def ocp1 begin
    tf ∈ R, variable
    t ∈ [ 0, tf ], time
    0.1 ≤ tf ≤ 10
    x ∈ R^3, state
    r = x_1
    v = x_2
    m = x_3
    x(0) == [1, 0, 1]
    m(tf) == 0.6
    v(t) <= 0.1
    u ∈ R, control
    0 ≤ u(t) ≤ 1
    Cd = 310
    β = 500
    Tmax = 3.5
    b = 2
    D = Cd * v^2 * exp(-β*(r - 1))
    F0: x -> [v, -D/m - 1/r^2, 0]
    F1: x -> [ 0, Tmax/m, -b*Tmax ]
    ẋ(t) == F0(x(t)) + u(t)*F1(x(t)) 
    tf → min
end
```

We can then solve the problem directly
```@example main
sol1 = solve(ocp1, print_level=5)
```
and plot the solution
```@example main
plot(sol1)
```

+++ variantes initial guess: constant, fonctionnel, mixed, warmstart

+++ variantes appel: transcription puis solve docp, avec init changee dans le docp (ex: mixed) et passee au solve (ex: warmstart)

