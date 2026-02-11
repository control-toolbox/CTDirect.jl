println("testing: nlp solver/modeler options")

if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
if !isdefined(Main, :goddard)
    include("../problems/goddard.jl")
end


# modeler / solver combinations
@testset verbose = true showtiming = true ":adnlp :exa :ipopt :madnlp" begin
    problem_list = [simple_integrator(), goddard2()]
    modeler_list = [:adnlp, :exa]
    solver_list = [:ipopt, :madnlp]
    for prob in problem_list
        for modeler in modeler_list
            for solver in solver_list
                test_problem(prob; modeler=modeler, solver=solver)
            end
        end
    end
end

# backends for ADNLPModels
@testset verbose = true showtiming = true ":adnlp_backends" begin
    prob = goddard()
    test_problem(prob)
    test_problem(prob; adnlp_backend=:default)
    test_problem(prob; adnlp_backend=:manual)
    test_problem(prob; scheme=:midpoint, adnlp_backend=:manual)
    test_problem(prob; scheme=:gauss_legendre_2, adnlp_backend=:manual)
end


# solution building
if !isdefined(Main, :double_integrator_minenergy)
    include("../problems/double_integrator.jl")
end
prob = double_integrator_minenergy(1)
prob2 = double_integrator_minenergy2(1)
x_opt = t -> [6 * (t^2 / 2 - t^3 / 3), 6 * (t - t^2)]
u_opt = t -> 6 - 12 * t
p_opt = t -> [24, 12 - 24 * t]

@testset verbose = true showtiming = true ":analytic_solution :ipopt :adnlp" begin
    sol = solve_problem(prob)
    T = time_grid(sol)
    @test isapprox(x_opt.(T), state(sol).(T), rtol=1e-2)
    @test isapprox(u_opt.(T), control(sol).(T), rtol=1e-2)
    @test isapprox(p_opt.(T), costate(sol).(T), rtol=1e-2)
end

@testset verbose = true showtiming = true ":analytic_solution :madnlp :adnlp" begin
    sol = solve_problem(prob; solver=:madnlp)
    T = time_grid(sol)
    @test isapprox(x_opt.(T), state(sol).(T), rtol=1e-2)
    @test isapprox(u_opt.(T), control(sol).(T), rtol=1e-2)
    @test isapprox(p_opt.(T), costate(sol).(T), rtol=1e-2)
end

@testset verbose = true showtiming = true ":analytic_solution :ipopt :exa" begin
    sol = solve_problem(prob2; modeler=:exa)
    T = time_grid(sol)
    @test isapprox(x_opt.(T), state(sol).(T), rtol=1e-2)
    @test isapprox(u_opt.(T[1:end-1]), control(sol).(T[1:end-1]), rtol=1e-2)
    @test isapprox(p_opt.(T), costate(sol).(T), rtol=1e-2)
end

@testset verbose = true showtiming = true ":analytic_solution :madnlp :exa" begin
    sol = solve_problem(prob2; solver=:madnlp, modeler=:exa)
    T = time_grid(sol)
    @test isapprox(x_opt.(T), state(sol).(T), rtol=1e-2)
    @test isapprox(u_opt.(T[1:end-1]), control(sol).(T[1:end-1]), rtol=1e-2)
    @test isapprox(p_opt.(T), costate(sol).(T), rtol=1e-2)
end


