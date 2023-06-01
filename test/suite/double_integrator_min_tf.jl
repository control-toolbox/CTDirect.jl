using CTDirect
using CTBase


# min tf
ocp1 = Model(variable=true)
state!(ocp1, 2)
control!(ocp1, 1)
variable!(ocp1, 1)
time!(ocp1, 0, Index(1))
constraint!(ocp1, :initial, [0,0], :initial_constraint)
constraint!(ocp1, :final, [1,0], :final_constraint)
constraint!(ocp1, :control, -1, 1, :control_constraint)
constraint!(ocp1, :variable, 0.1, 10, :variable_constraint)
dynamics!(ocp1, (x, u, v) ->  [x[2], u])
objective!(ocp1, :mayer, (x0, xf, v) -> v) # scalar variable
#objective!(ocp3, :mayer, (x0, xf, v) -> - v[1]) # max t0 free to and tf ok
#objective!(ocp3, :mayer, (x0, xf, v) -> v[1], :max) # max t0 free to and tf ok
sol1 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
@testset verbose = true showtiming = true ":double_integrator :min_tf" begin
    @test sol1.objective ≈ 2.0 rtol=1e-2
end


# min tf (lagrange)
ocp2 = Model(variable=true)
state!(ocp2, 2)
control!(ocp2, 1)
variable!(ocp2, 1)
time!(ocp2, 0, Index(1))
constraint!(ocp2, :initial, [0,0], :initial_constraint)
constraint!(ocp2, :final, [1,0], :final_constraint)
constraint!(ocp2, :control, -1, 1, :control_constraint)
constraint!(ocp2, :variable, 0.1, 10, :variable_constraint)
dynamics!(ocp2, (x, u, v) ->  [x[2], u])
objective!(ocp2, :lagrange, (x, u, v) -> 1)
sol2 = solve(ocp2, grid_size=100, print_level=0, tol=1e-12)
@testset verbose = true showtiming = true ":double_integrator :min_tf :lagrange" begin
    @test sol2.objective ≈ 2.0 rtol=1e-2
end


# min tf (vectorial)
ocp3 = Model(variable=true)
state!(ocp3, 2)
control!(ocp3, 1)
variable!(ocp3, 1)
time!(ocp3, 0, Index(1))
constraint!(ocp3, :initial, [0,0], :initial_constraint)
constraint!(ocp3, :final, [1,0], :final_constraint)
constraint!(ocp3, :control, [-1], [1], :control_constraint)
constraint!(ocp3, :variable, [0.1], [10], :variable_constraint)
dynamics!(ocp3, (x, u, v) ->  [x[2], u[1]])
objective!(ocp3, :mayer, (x0, xf, v) -> v[1])
sol3 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
@testset verbose = true showtiming = true ":double_integrator :min_tf :vectorial" begin
    @test sol3.objective ≈ 2.0 rtol=1e-2
end


# max t0 (free t0 and tf)
ocp4 = Model(variable=true)
state!(ocp4, 2)
control!(ocp4, 1)
variable!(ocp4, 2)
time!(ocp4, Index(1), Index(2))
constraint!(ocp4, :initial, [0,0], :initial_constraint)
constraint!(ocp4, :final, [1,0], :final_constraint)
constraint!(ocp4, :control, -1, 1, :control_constraint)
constraint!(ocp4, :variable, [0.1, 0.1], [10, 10], :variable_constraint)
dynamics!(ocp4, (x, u, v) ->  [x[2], u])
objective!(ocp4, :mayer, (x0, xf, v) -> v[1], :max)
sol4 = solve(ocp4, grid_size=100, print_level=0, tol=1e-12)
@testset verbose = true showtiming = true ":double_integrator :max_t0" begin
    @test sol4.objective ≈ 8.0 rtol=1e-2
end