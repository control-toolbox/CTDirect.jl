println("Test: misc")

@testset verbose = true ":default :direct" begin
    @test CTDirect.__grid_size() isa Integer
    @test isnothing(CTDirect.__time_grid())
    @test CTDirect.__tolerance() < 1
    @test CTDirect.__max_iterations() isa Integer
end

@testset verbose = true ":default :ipopt" begin
    @test CTDirect.__ipopt_print_level() isa Integer
    @test CTDirect.__ipopt_print_level() ≤ 12
    @test CTDirect.__ipopt_print_level() ≥ 0
    @test CTDirect.__ipopt_mu_strategy() isa String
    @test CTDirect.__ipopt_linear_solver() isa String
end

# simple integrator min energy dual control
if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
ocp = simple_integrator().ocp
sol0 = direct_solve(ocp, display = false)

# test save / load solution in JLD2 format
@testset verbose = true showtiming = true ":save_load :JLD2" begin
    save(sol0, filename_prefix = "solution_test")
    sol_reloaded = load("solution_test")
    @test sol0.objective == sol_reloaded.objective
end

# test export / read solution in JSON format
@testset verbose = true showtiming = true ":export_read :JSON" begin
    export_ocp_solution(sol0, filename_prefix = "solution_test")
    sol_reloaded = import_ocp_solution("solution_test")
    @test sol0.objective == sol_reloaded.objective
end
