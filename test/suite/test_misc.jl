println("Test: misc")

@testset verbose = true ":default_direct" begin
    @test CTDirect.__grid_size() isa Integer
    @test isnothing(CTDirect.__time_grid())
    @test CTDirect.__tolerance() < 1
    @test CTDirect.__max_iterations() isa Integer
end

@testset verbose = true ":default_ipopt" begin
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
