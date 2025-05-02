println("testing: nlp options")

if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
prob = simple_integrator()
ocp = prob.ocp
obj = prob.obj

@testset verbose = true showtiming = true ":methods" begin
    @test CTDirect.is_solvable(ocp)
    @test (:adnlp, :ipopt) in CTDirect.available_methods()
    @test (:adnlp, :madnlp) in CTDirect.available_methods()
end

# AD backends
@testset verbose = true showtiming = true ":AD_backends" begin
    sol = solve(ocp, display=false)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = solve(ocp, display=false, adnlp_backend=:default)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = solve(ocp, display=false, adnlp_backend=:manual)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = solve(ocp, display=false, disc_method=:midpoint, adnlp_backend=:manual)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = solve(ocp, display=false, disc_method=:gauss_legendre_2, adnlp_backend=:manual)
    @test objective(sol) ≈ obj rtol = 1e-2
end


# DOCP solving
@testset verbose = true showtiming = true ":solve_docp" begin
    docp, nlp = direct_transcription(ocp)
    solver_backend = CTDirect.IpoptBackend()
    dsol = CTDirect.solve_docp(solver_backend, docp, nlp, display=false)
    sol = CTDirect.build_OCP_solution(docp, dsol)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = CTDirect.build_OCP_solution(docp, primal=dsol.solution)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = CTDirect.build_OCP_solution(docp, primal=dsol.solution, dual=dsol.multipliers)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = CTDirect.build_OCP_solution(docp, primal=dsol.solution, dual=dsol.multipliers, mult_LB=dsol.multipliers_L, mult_UB=dsol.multipliers_U)
    @test objective(sol) ≈ obj rtol = 1e-2
end


@testset verbose = true showtiming = true ":solve_docp :madnlp :gl2" begin
    docp, nlp = direct_transcription(ocp, disc_method=:gauss_legendre_2)
    solver_backend = CTDirect.MadNLPBackend()
    dsol = CTDirect.solve_docp(solver_backend, docp, nlp, display=false)
    sol = CTDirect.build_OCP_solution(docp, dsol)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = CTDirect.build_OCP_solution(docp, primal=dsol.solution)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = CTDirect.build_OCP_solution(docp, primal=dsol.solution, dual=dsol.multipliers)
    @test objective(sol) ≈ obj rtol = 1e-2
    sol = CTDirect.build_OCP_solution(docp, primal=dsol.solution, dual=dsol.multipliers, mult_LB=dsol.multipliers_L, mult_UB=dsol.multipliers_U)
    @test objective(sol) ≈ obj rtol = 1e-2
end

# solution building
if !isdefined(Main, :double_integrator_minenergy)
    include("../problems/double_integrator.jl")
end
ocp = double_integrator_minenergy(1).ocp
x_opt = t -> [6 * (t^2 / 2 - t^3 / 3), 6 * (t - t^2)]
u_opt = t -> 6 - 12 * t
p_opt = t -> [24, 12 - 24 * t]

@testset verbose = true showtiming = true ":analytic_solution :ipopt" begin
    sol = solve(ocp, display=false)
    T = time_grid(sol)
    @test isapprox(x_opt.(T), state(sol).(T), rtol=1e-2)
    @test isapprox(u_opt.(T), control(sol).(T), rtol=1e-2)
    @test isapprox(p_opt.(T), costate(sol).(T), rtol=1e-2)
end

@testset verbose = true showtiming = true ":analytic_solution :madnlp" begin
    sol = solve(ocp, :madnlp, display=false)
    T = time_grid(sol)
    @test isapprox(x_opt.(T), state(sol).(T), rtol=1e-2)
    @test isapprox(u_opt.(T), control(sol).(T), rtol=1e-2)
    @test isapprox(p_opt.(T), costate(sol).(T), rtol=1e-2)
end


# setting the initial guess at the DOCP level
prob = double_integrator_mintf()
ocp = prob.ocp
sol0 = solve(ocp, display=false)
docp, nlp = direct_transcription(ocp)
solver_backend = CTDirect.IpoptBackend()
v_const = 0.15
t_vec = [0, 0.1, v_const]
x_vec = [[0, 0], [1, 2], [5, -1]]
u_func = t -> (cos(10 * t) + 1) * 0.5
# mixed init
@testset verbose = true showtiming = true ":docp_mixed_init" begin
    set_initial_guess(
        docp,
        nlp,
        (time=t_vec, state=x_vec, control=u_func, variable=v_const),
    )
    dsol = CTDirect.solve_docp(solver_backend, docp, nlp, display=false, max_iter=maxiter)
    sol = CTDirect.build_OCP_solution(docp, dsol)
    T = time_grid(sol)
    @test isapprox(state(sol).(t_vec), x_vec, rtol=1e-2)
    @test isapprox(control(sol).(T), u_func.(T), rtol=1e-2)
    @test variable(sol) == v_const
end
# warm start
@testset verbose = true showtiming = true ":docp_warm_start" begin
    set_initial_guess(docp, nlp, sol0)
    dsol = CTDirect.solve_docp(solver_backend, docp, nlp, display=false, max_iter=maxiter)
    sol = CTDirect.build_OCP_solution(docp, dsol)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), state(sol0).(T), rtol=1e-2)
    @test isapprox(control(sol).(T), control(sol0).(T), rtol=1e-2)
    @test variable(sol) == variable(sol0)
end
