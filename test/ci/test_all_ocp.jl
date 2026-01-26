println("testing: OCP definitions")

# beam
if !isdefined(Main, :beam)
    include("../problems/beam.jl")
end
@testset verbose = true showtiming = true ":beam" begin
    test_problem(beam())
end

# bolza, non-autonomous mayer term, tf in dynamics
if !isdefined(Main, :bolza_freetf)
    include("../problems/bolza.jl")
end
@testset verbose = true showtiming = true ":bolza :tf_in_dyn_and_cost" begin
    test_problem(bolza_freetf())
end

# double integrator min tf / max t0 (free t0 and tf)
if !isdefined(Main, :double_integrator_mintf)
    include("../problems/double_integrator.jl")
end
@testset verbose = true showtiming = true ":double_integrator" begin
    test_problem(double_integrator_mintf())
    test_problem(double_integrator_freet0tf())
    test_problem(double_integrator_nobounds())
end

# electric vehicle
if !isdefined(Main, :electric_vehicle)
    include("../problems/electric_vehicle.jl")
end
@testset verbose = true showtiming = true ":electric_vehicle" begin
    test_problem(electric_vehicle())
end

# fuller
if !isdefined(Main, :fuller)
    include("../problems/fuller.jl")
end
@testset verbose = true showtiming = true ":fuller" begin
    test_problem(fuller())
end

# glider
if !isdefined(Main, :glider)
    include("../problems/glider.jl")
end
@testset verbose = true showtiming = true ":glider" begin
    test_problem(glider())
end

# goddard max rf (with all constraints version)
if !isdefined(Main, :goddard)
    include("../problems/goddard.jl")
end
@testset verbose = true showtiming = true ":goddard :max_rf" begin
    test_problem(goddard())
    test_problem(goddard_all())
end

# insurance (nb. requires final control for CV, mixed constraints)
if !isdefined(Main, :insurance)
    include("../problems/insurance.jl")
end
@testset verbose = true showtiming = true ":insurance" begin
    test_problem(insurance(), scheme=:trapeze)
end

# jackson
if !isdefined(Main, :jackson)
    include("../problems/jackson.jl")
end
@testset verbose = true showtiming = true ":jackson" begin
    test_problem(jackson())
end

# moonlander
if !isdefined(Main, :moonlander)
    include("../problems/moonlander.jl")
end
@testset verbose = true showtiming = true ":moonlander" begin
    test_problem(moonlander(), adnlp_backend=:manual)
end

# quadrotor
if !isdefined(Main, :quadrotor)
    include("../problems/quadrotor.jl")
end
@testset verbose = true showtiming = true ":quadrotor" begin
    test_problem(moonlander(), adnlp_backend=:manual)
end

# robbins
if !isdefined(Main, :robbins)
    include("../problems/robbins.jl")
end
@testset verbose = true showtiming = true ":robbins" begin
    test_problem(robbins())
end

# simple integrator
if !isdefined(Main, :simple_integrator)
    include("../problems/simple_integrator.jl")
end
@testset verbose = true showtiming = true ":simple_integrator" begin
    test_problem(simple_integrator())
end

# space shuttle
if !isdefined(Main, :space_shuttle)
    include("../problems/space_shuttle.jl")
end
@testset verbose = true showtiming = true ":space_shuttle" begin
    test_problem(space_shuttle())
end

# truck trailer
if !isdefined(Main, :truck_trailer)
    include("../problems/truck_trailer.jl")
end
@testset verbose = true showtiming = true ":truck_trailer" begin
    test_problem(truck_trailer())
end

# vanderpol
if !isdefined(Main, :vanderpol)
    include("../problems/vanderpol.jl")
end
@testset verbose = true showtiming = true ":vanderpol" begin
    test_problem(vanderpol())
end
