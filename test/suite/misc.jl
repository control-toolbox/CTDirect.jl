println("Test: misc")

# simple integrator min energy
# control split as positive/negative parts for m=2 tets case
ocp = Model()
state!(ocp, 1)
control!(ocp, 2)
time!(ocp, t0=0, tf=1)
constraint!(ocp, :initial, lb=-1, ub=-1)
constraint!(ocp, :final, lb=0, ub=0)
constraint!(ocp, :control, lb=[0,0], ub=[Inf, Inf])
dynamics!(ocp, (x, u) -> -x - u[1] + u[2])
objective!(ocp, :lagrange, (x, u) -> (u[1]+u[2])^2)
sol0 = solve(ocp, print_level=0)

@testset verbose = true showtiming = true ":control_dim_2" begin
    @test (:adnlp,:ipopt) in available_methods()
    @test is_solvable(ocp)
end

@testset verbose = true showtiming = true ":docp_solve" begin
    docp = directTranscription(ocp)
    dsol = solve(docp, print_level=0, tol=1e-12)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.objective ≈ 0.313 rtol=1e-2
end

@testset verbose = true showtiming = true ":docp_solve :warm_start" begin
    docp = directTranscription(ocp, init=sol0)
    dsol = solve(docp, print_level=0, tol=1e-12)
    sol = OCPSolutionFromDOCP(docp, dsol)
    @test sol.iterations == 5
end

# test NLP getter
docp = directTranscription(ocp)
nlp = getNLP(docp)

# test OCPSolutionFromDOCP_raw
#dsol = solve(docp, print_level=0, tol=1e-12)
#sol_raw = OCPSolutionFromDOCP_raw(docp, dsol.solution)

# test save / load solution in JLD2 format
@testset verbose = true showtiming = true ":save_load :JLD2" begin
    save_OCP_solution(sol0, filename_prefix="solution_test")
    sol_reloaded = load_OCP_solution("solution_test")
    @test sol0.objective == sol_reloaded.objective
end

# test export / read solution in JSON format
@testset verbose = true showtiming = true ":export_read :JSON" begin
    export_OCP_solution(sol0, filename_prefix="solution_test")
    sol_reloaded = read_OCP_solution("solution_test")
    @test sol0.objective == sol_reloaded.objective
end

