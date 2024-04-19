# Discrete continuation

Using the warm start option, it is easy to implement a basic discrete continuation method, where a sequence of problems is solved using each solution as initial guess for the next problem.
This usually gives better and faster convergence than solving each problem with the same initial guess.
We illustrate this on the Goddard problem already presented in the tutorial.

```@setup main
using Printf
using Statistics
using Plots
```

```@example main
using CTDirect
using CTBase

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

## Continuation on speed constraint

Let us solve a sequence of problems with an increasingly strict constraint on the maximal speed.
```@example main
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

We now plot the objective with respect to the speed limit, as well as a comparison of the solutions for the unconstrained case and the vmax=0.05 case.

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

## Continuation on maximal thrust

As a second example, we now solve the problem for a decreasing maximal thrust Tmax. first we reset the speed constraint.

```@example main
remove_constraint!(ocp, :speed_limit)
vmax = 0.1
constraint!(ocp, :state, Index(2), 0, vmax, :speed_limit)
sol = solve(ocp, print_level=0)
sol0 = sol
nothing # hide
```

Then we perform the continuation

```@example main
Tmax_list = []
obj_list = []
print("Tmax ")
for Tmax_local=3.5:-0.5:1
    global Tmax = Tmax_local 
    print(Tmax," ")   
    global sol = solve(ocp, print_level=0, init=OptimalControlInit(sol))
    push!(Tmax_list, Tmax)
    push!(obj_list, sol.objective)
end 
```

And plot the results as before

```@example main
pobj = plot(Tmax_list, obj_list, label="r(tf)",seriestype=:scatter)
xlabel!("Maximal thrust (Tmax)")
ylabel!("Maximal altitude r(tf)")
plot(sol0)
p = plot!(sol)
plot(pobj, p, layout=2)
```
