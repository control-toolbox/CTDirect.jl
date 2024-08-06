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

# +++ move this in a new test_nlp.jl
@testset verbose = true showtiming = true ":control_dim_2" begin
    @test (:adnlp,:ipopt) in available_methods()
    @test (:adnlp,:madnlp) in available_methods()
    @test is_solvable(ocp)
end

@testset verbose = true showtiming = true ":docp_solve" begin
    docp, nlp = direct_transcription(ocp)
    solver = IpoptSolver(nlp)
    dsol = CTDirect.solve_docp(solver, docp, nlp, print_level=0, tol=1e-12)
    sol = build_solution(docp, dsol)
    @test sol.objective ≈ 0.313 rtol=1e-2
end

@testset verbose = true showtiming = true ":docp_solve :warm_start" begin
    docp, nlp = direct_transcription(ocp, init=sol0)
    solver = IpoptSolver(nlp)
    dsol = CTDirect.solve_docp(solver, docp, nlp, print_level=0, tol=1e-12)
    sol = build_solution(docp, dsol)
    @test sol.iterations == 5
end

# build solution from NLP (+++ add more tests ?)
@testset verbose = true showtiming = true ":build_solution_from_nlp" begin
    docp, nlp = direct_transcription(ocp, init=sol0)
    solver = IpoptSolver(nlp)
    dsol = CTDirect.solve_docp(solver, docp, nlp, print_level=0, tol=1e-12)
    sol = build_solution(docp, primal=dsol.solution)
    @test sol.objective ≈ 0.313 rtol=1e-2
    sol = build_solution(docp, primal=dsol.solution, dual=dsol.multipliers)
    @test sol.objective ≈ 0.313 rtol=1e-2
end

# test save / load solution in JLD2 format
@testset verbose = true showtiming = true ":save_load :JLD2" begin
    save(sol0, filename_prefix="solution_test")
    sol_reloaded = load("solution_test")
    @test sol0.objective == sol_reloaded.objective
end

# test export / read solution in JSON format
@testset verbose = true showtiming = true ":export_read :JSON" begin
    export_ocp_solution(sol0, filename_prefix="solution_test")
    sol_reloaded = import_ocp_solution("solution_test")
    @test sol0.objective == sol_reloaded.objective
end
