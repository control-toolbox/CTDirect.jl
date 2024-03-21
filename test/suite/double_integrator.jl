using CTDirect

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
#sol1 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
docp1 = DirectTranscription(ocp1, grid_size=100);
sol1 = solveDOCP(docp1, print_level=0, tol=1e-12);
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
#sol2 = solve(ocp2, grid_size=100, print_level=0, tol=1e-12)
docp2 = DirectTranscription(ocp2, grid_size=100);
sol2 = solveDOCP(docp2, print_level=0, tol=1e-12);
@testset verbose = true showtiming = true ":double_integrator :min_tf :lagrange" begin
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
#sol3 = solve(ocp1, grid_size=100, print_level=0, tol=1e-12)
docp3 = DirectTranscription(ocp3, grid_size=100);
sol3 = solveDOCP(docp3, print_level=0, tol=1e-12);
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
#sol4 = solve(ocp4, grid_size=100, print_level=0, tol=1e-12)
docp4 = DirectTranscription(ocp4, grid_size=100);
sol4 = solveDOCP(docp4, print_level=0, tol=1e-12);
@testset verbose = true showtiming = true ":double_integrator :max_t0" begin
    @test sol4.objective ≈ 8.0 rtol=1e-2
end


# min tf, abstract definition
@def ocp begin
    tf ∈ R, variable
    t ∈ [ 0, tf ], time
    x ∈ R², state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    x(0) == [ 0, 0 ]
    x(tf) == [ 1, 0 ]
    0.1 ≤ tf ≤ Inf 
    ẋ(t) == [ x₂(t), u(t) ] 
    tf → min
end
#sol = solve(ocp, grid_size=100, print_level=0, tol=1e-12)
docp = DirectTranscription(ocp, grid_size=100);
sol = solveDOCP(docp, print_level=0, tol=1e-12);
@testset verbose = true showtiming = true ":double_integrator :min_tf :abstract" begin
    @test is_solvable(ocp)
    @test sol.objective ≈ 2.0 rtol=1e-2
end


# min energy dual control
ocp5 = Model()
state!(ocp5, 2)
control!(ocp5, 2)
time!(ocp5, 0, 5)
constraint!(ocp5, :initial, [0,0], :initial_constraint)
constraint!(ocp5, :final, [1,0], :final_constraint)
constraint!(ocp5, :control, 1:2, [0,0], [1,1], :control_box) # [u_, u+]
dynamics!(ocp5, (x, u) ->  [x[2], -u[1] + u[2]])
objective!(ocp5, :lagrange, (x, u) -> u[1]*u[1] + u[2]*u[2])
#sol5 = solve(ocp5, grid_size=50, print_level=0, tol=1e-12)
docp5 = DirectTranscription(ocp5, grid_size=100);
sol5 = solveDOCP(docp5, print_level=0, tol=1e-12);
@testset verbose = true showtiming = true ":double_integrator :min_energy" begin
    @test sol5.objective ≈ 9.6e-2 rtol=1e-2
end