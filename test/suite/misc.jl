using CTDirect

println("Test: misc")

# simple integrator min energy
# control split as positive/negative parts for m=2 tets case
ocp = Model()
state!(ocp, 1)
control!(ocp, 2)
time!(ocp, [0, 1])
constraint!(ocp, :initial, -1, :initial_constraint)
constraint!(ocp, :final, 0, :final_constraint)
constraint!(ocp, :control, 1:2, [0,0], [Inf, Inf], :positive_controls)
dynamics!(ocp, (x, u) -> -x - u[1] + u[2])
objective!(ocp, :lagrange, (x, u) -> (u[1]+u[2])^2)
sol0 = solve(ocp, print_level=0)

@testset verbose = true showtiming = true ":control_dim_2" begin
    @test (:adnlp,:ipopt) in available_methods()
    @test is_solvable(ocp)
end

@testset verbose = true showtiming = true ":docp_solve" begin
    docp = directTranscription(ocp, grid_size=100)
    sol, dsol = solve(docp, print_level=0, tol=1e-12)
    @test sol.objective ≈ 0.313 rtol=1e-2
end

@testset verbose = true showtiming = true ":docp_solve :warm_start" begin
    docp2 = directTranscription(ocp, grid_size=100, init=sol0)
    sol2, dsol2 = solve(docp2, print_level=0, tol=1e-12)
    @test sol2.iterations == 5
end

# test NLP getter
docp = directTranscription(ocp, grid_size=100)
nlp = getNLP(docp)

# test OCPSolutionFromDOCP_raw
sol, dsol = solve(docp, print_level=0, tol=1e-12)
sol_raw = OCPSolutionFromDOCP_raw(docp, dsol.solution)

