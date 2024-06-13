# tests some objective options, variable tf
println("Test: double integrator")

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
objective!(ocp1, :mayer, (x0, xf, v) -> v)

@testset verbose = true showtiming = true ":min_tf :mayer" begin
    sol1 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
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

@testset verbose = true showtiming = true ":min_tf :lagrange" begin
    sol2 = solve(ocp2, grid_size=100, print_level=0, tol=1e-12)
    @test sol2.objective ≈ 2.0 rtol=1e-2
end


# min tf (vector)
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

@testset verbose = true showtiming = true ":min_tf :mayer :vector" begin
    sol3 = solve(ocp3, grid_size=100, print_level=0, tol=1e-12)
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
constraint!(ocp4, :variable, v -> v[2]-v[1], 0.1, Inf)
dynamics!(ocp4, (x, u, v) ->  [x[2], u])
objective!(ocp4, :mayer, (x0, xf, v) -> v[1], :max)

@testset verbose = true showtiming = true ":max_t0" begin
    sol4 = solve(ocp4, grid_size=100, print_level=0, tol=1e-12)
    @test sol4.objective ≈ 8.0 rtol=1e-2
end

@testset verbose = true showtiming = true ":max_t0 :non_uniform_grid" begin
    sol5 = solve(ocp4, time_grid=[0,0.1,0.6,0.95,1], print_level=0, tol=1e-12)
    @test sol5.objective ≈ 7.48 rtol=1e-2
end
