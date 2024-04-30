# Discrete continuation

Using the warm start option, it is easy to implement a basic discrete continuation method, where a sequence of problems is solved using each solution as initial guess for the next problem.
This usually gives better and faster convergence than solving each problem with the same initial guess, and is a way to handle problems that require a good initial guess.


## Continuation on parametric OCP

The most compact syntax to perform a discrete continuation is to use a function that returns the OCP for a given value of the continuation parameter, and solve a sequence of these problems. We illustrate this on a very basic double integrator with increasing fixed final time.

```@setup main
using Printf
using Statistics
using Plots
```

First we write a function that returns the OCP for a given final time.

```@example main
using CTDirect
using CTBase

function ocp_T(T)
    @def ocp begin
        t ∈ [ 0, T ], time
        x ∈ R², state
        u ∈ R, control
        q = x₁
        v = x₂
        q(0) == 0
        v(0) == 0
        q(T) == 1
        v(T) == 0
        ẋ(t) == [ v(t), u(t) ]
        ∫(u(t)^2) → min
    end
    return ocp
end
nothing # hide
```

Then we perform the continuation with a simple *for* loop, using each solution to initialize the next problem.

```@example main
init1 = OptimalControlInit()
iter_list = []
for T=1:5
    ocp1 = ocp_T(T) 
    sol1 = solve(ocp1, print_level=0, init=init1)
    global init1 = OptimalControlInit(sol1)
    @printf("T %.2f objective %.6f iterations %d\n", T, sol1.objective, sol1.iterations)
    push!(iter_list, sol1.iterations)
end
```

## Continuation on global variable

As a second example, we show how to avoid redefining a new OCP each time, and modify the original one instead.
More precisely we now solve a Goddard problem for a decreasing maximal thrust. If we store the value for *Tmax* in a global variable, we can simply modify this variable and keep the same OCP problem during the continuation.

Let us first define the Goddard problem (note that the formulation below illustrates all the possible constraints types, and the problem could be defined in a more compact way).

```@example main
Cd = 310
Tmax = 3.5
β = 500
b = 2
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end

ocp = Model(variable=true)
state!(ocp, 3)
control!(ocp, 1)
variable!(ocp, 1)
time!(ocp, 0, Index(1))
constraint!(ocp, :initial, [1,0,1], :initial_constraint)
constraint!(ocp, :final, Index(3), 0.6, :final_constraint)
constraint!(ocp, :state, 1:2:3, [1,0.6], [1.2,1], :state_box)
constraint!(ocp, :control, Index(1), 0, 1, :control_box)
constraint!(ocp, :variable, Index(1), 0.01, Inf, :variable_box)
constraint!(ocp, :state, Index(2), 0, Inf, :speed_limit)
objective!(ocp, :mayer, (x0, xf, v) -> xf[1], :max)
dynamics!(ocp, (x, u, v) -> F0(x) + u*F1(x) )

sol0 = solve(ocp, print_level=0)
sol = sol0
@printf("Objective for reference solution %.6f\n", sol0.objective)
```

Then we perform the continuation on the maximal thrust.

```@example main
Tmax_list = []
obj_list = []
for Tmax_local=3.5:-0.5:1
    global Tmax = Tmax_local  
    global sol = solve(ocp, print_level=0, init=OptimalControlInit(sol))
    @printf("Tmax %.2f objective %.6f iterations %d\n", Tmax, sol.objective, sol.iterations)
    push!(Tmax_list, Tmax)
    push!(obj_list, sol.objective)
end 
```

We plot now the objective w.r.t the maximal thrust, as well as both solutions for *Tmax*=3.5 and *Tmax*=1.

```@example main
pobj = plot(Tmax_list, obj_list, label="r(tf)",seriestype=:scatter)
xlabel!("Maximal thrust (Tmax)")
ylabel!("Maximal altitude r(tf)")
plot(sol0)
p = plot!(sol)
plot(pobj, p, layout=2)
```


## Manual constraint redefinition

Here we illustrate a slightly more involved way of modifying the OCP problem during the continuation.
Instead of just updating a global variable as before, we now remove and redefine one of the constraints (maximal speed). 
```@example main
global Tmax = 3.5
vmax_list = []
obj_list = []
iter_list = []
print("vmax ")
for vmax=0.15:-0.01:0.05
    print(vmax," ")
    remove_constraint!(ocp, :speed_limit)
    constraint!(ocp, :state, Index(2), 0, vmax, :speed_limit)
    global sol = solve(ocp, print_level=0, init=OptimalControlInit(sol))
    push!(vmax_list, vmax)
    push!(obj_list, sol.objective)
    push!(iter_list, sol.iterations)
end
@printf("\nAverage iterations %d\n", mean(iter_list))
```

We now plot the objective with respect to the speed limit, as well as a comparison of the solutions for the unconstrained case and the *vmax*=0.05 case.

```@example main
pobj = plot(vmax_list, obj_list, label="r(tf)",seriestype=:scatter)
xlabel!("Speed limit (vmax)")
ylabel!("Maximal altitude r(tf)")
plot(sol0)
p = plot!(sol)
plot(pobj, p, layout=2)
```

We can compare with solving each problem with the default initial guess, which here gives the same solutions but takes more iterations overall.

```@example main
iter_list = []
for vmax=0.15:-0.01:0.05
    print(vmax," ")
    remove_constraint!(ocp, :speed_limit)
    constraint!(ocp, :state, Index(2), 0, vmax, :speed_limit)
    global sol = solve(ocp, print_level=0)   
    push!(iter_list, sol.iterations)
end
@printf("\nAverage iterations %d\n", mean(iter_list))
```