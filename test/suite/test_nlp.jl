println("testing: nlp options")

if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
prob = simple_integrator()
ocp = prob.ocp

@testset verbose = true showtiming = true ":control_dim_2" begin
    @test is_solvable(ocp)
    @test (:adnlp, :ipopt) in available_methods()
    @test (:adnlp, :madnlp) in available_methods()
end

@testset verbose = true showtiming = true ":solve_docp" begin
    docp, nlp = direct_transcription(ocp)
    solver_backend = CTDirect.IpoptBackend()
    dsol = CTDirect.solve_docp(solver_backend, docp, nlp, display=false)
    sol = OptimalControlSolution(docp, dsol)
    @test sol.objective ≈ prob.obj rtol = 1e-2
    sol = OptimalControlSolution(docp, primal=dsol.solution)
    @test sol.objective ≈ prob.obj rtol = 1e-2
    sol = OptimalControlSolution(docp, primal=dsol.solution, dual=dsol.multipliers)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

@testset verbose = true showtiming = true ":solve_docp :midpoint" begin
    docp, nlp = direct_transcription(ocp, disc_method=:midpoint)
    solver_backend = CTDirect.IpoptBackend()
    dsol = CTDirect.solve_docp(solver_backend, docp, nlp, display=false)
    sol = OptimalControlSolution(docp, dsol)
    @test sol.objective ≈ prob.obj rtol = 1e-2
    sol = OptimalControlSolution(docp, primal=dsol.solution)
    @test sol.objective ≈ prob.obj rtol = 1e-2
    sol = OptimalControlSolution(docp, primal=dsol.solution, dual=dsol.multipliers)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

@testset verbose = true showtiming = true ":solve_docp :madnlp" begin
    docp, nlp = direct_transcription(ocp)
    solver_backend = CTDirect.MadNLPBackend()
    dsol = CTDirect.solve_docp(solver_backend, docp, nlp, display=false)
    sol = OptimalControlSolution(docp, dsol)
    @test sol.objective ≈ prob.obj rtol = 1e-2
    sol = OptimalControlSolution(docp, primal=dsol.solution)
    @test sol.objective ≈ prob.obj rtol = 1e-2
    sol = OptimalControlSolution(docp, primal=dsol.solution, dual=dsol.multipliers)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# check solution building
if !isdefined(Main, :double_integrator_minenergy)
    include("../problems/double_integrator.jl")
end
ocp = double_integrator_minenergy(1).ocp
x_opt = t -> [6 * (t^2 / 2 - t^3 / 3), 6 * (t - t^2)]
u_opt = t -> 6 - 12 * t
p_opt = t -> [24, 12 - 24 * t]

@testset verbose = true showtiming = true ":analytic_solution :ipopt" begin
    sol = direct_solve(ocp, display=false)
    T = sol.time_grid
    @test isapprox(x_opt.(T), sol.state.(T), rtol=1e-2)
    @test isapprox(u_opt.(T), sol.control.(T), rtol=1e-2)
    @test isapprox(p_opt.(T), sol.costate.(T), rtol=1e-2)
end

@testset verbose = true showtiming = true ":analytic_solution :madnlp" begin
    sol = direct_solve(ocp, :madnlp, display=false)
    T = sol.time_grid
    @test isapprox(x_opt.(T), sol.state.(T), rtol=1e-2)
    @test isapprox(u_opt.(T), sol.control.(T), rtol=1e-2)
    @test isapprox(p_opt.(T), sol.costate.(T), rtol=1e-2)
end
