println("testing: OCP definitions")

# beam
if !isdefined(Main, :beam)
    include("../problems/beam.jl")
end
@testset verbose = true showtiming = true ":beam" begin
    check_problem(beam(), display=false)
end

# bolza, non-autonomous mayer term, tf in dynamics
if !isdefined(Main, :bolza_freetf)
    include("../problems/bolza.jl")
end
@testset verbose = true showtiming = true ":bolza :tf_in_dyn_and_cost" begin
    check_problem(bolza_freetf(), display=false)
end

# double integrator min tf / max t0 (free t0 and tf)
if !isdefined(Main, :double_integrator_mintf)
    include("../problems/double_integrator.jl")
end
@testset verbose = true showtiming = true ":double_integrator :min_tf" begin
    check_problem(double_integrator_mintf(), display=false)
    check_problem(double_integrator_freet0tf(), display=false)
end

# fuller
if !isdefined(Main, :fuller)
    include("../problems/fuller.jl")
end
@testset verbose = true showtiming = true ":fuller" begin
    check_problem(fuller() , display=false)
end

# glider
if !isdefined(Main, :glider)
    include("../problems/glider.jl")
end
@testset verbose = true showtiming = true ":glider" begin
    check_problem(glider(), display=false)
end

# goddard max rf (with all constraints version)
if !isdefined(Main, :goddard)
    include("../problems/goddard.jl")
end
@testset verbose = true showtiming = true ":goddard :max_rf" begin
    check_problem(goddard(), display=false)
    check_problem(goddard_all() , display=false)
end

# insurance (nb. requires final control for CV, mixed constraints)
if !isdefined(Main, :insurance)
    include("../problems/insurance.jl")
end
@testset verbose = true showtiming = true ":insurance" begin
    check_problem(insurance() , display=false, disc_method=:trapeze)
end

# jackson
if !isdefined(Main, :jackson)
    include("../problems/jackson.jl")
end
@testset verbose = true showtiming = true ":jackson" begin
    check_problem(jackson(), display=false)
end

# moonlander
if !isdefined(Main, :moonlander)
    include("../problems/moonlander.jl")
end
@testset verbose = true showtiming = true ":moonlander" begin
    check_problem(moonlander(), display=false, adnlp_backend=:manual)
end

#= quadrotor
if !isdefined(Main, :quadrotor)
    include("../problems/quadrotor.jl")
end
@testset verbose = true showtiming = true ":quadrotor" begin
    prob = moonlander()
    sol = solve(prob.ocp, display=false, disc_method=:midpoint, adnlp_backend=:manual)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end=#

# robbins
if !isdefined(Main, :robbins)
    include("../problems/robbins.jl")
end
@testset verbose = true showtiming = true ":robbins" begin
    check_problem(robbins(), display=false)
end

# simple integrator
if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
@testset verbose = true showtiming = true ":simple_integrator" begin
    check_problem(simple_integrator(), display=false)
end

# space shuttle

# truck trailer

# vanderpol
if !isdefined(Main, :vanderpol)
    include("../problems/vanderpol.jl")
end
@testset verbose = true showtiming = true ":vanderpol" begin
    check_problem(vanderpol(), display=false)
end
