# runtests.jl
using Test
include("test_common.jl")

function test_problem(prob; kwargs...)

    sol = solve_problem(prob; kwargs...)
    @test CTModels.successful(sol)
    @test CTModels.iterations(sol) <= max_iter
    @test CTModels.constraints_violation(sol) <= tol
    @test sol.objective ≈ prob.obj rtol = 1e-2

end

# check local test suite
macro ignore(e)
    :()
end

# # run either usual test suite on CPU, or GPU tests only 
# @testset verbose = true showtiming = true "Test CTDirect" begin
#     if "GPU" in ARGS
#         # ExaModels tests only (GPU on moonshot workflow)
#         include("./suite/test_exa.jl")
#     else
#         # CPU: run all scripts in subfolder suite/
#         include.(filter(contains(r".jl$"), readdir("./suite"; join=true)))
#     end
# end


# new ci tests
const VERBOSE = true
const SHOWTIMING = true
@testset verbose = VERBOSE showtiming = SHOWTIMING "New tests for CTDirect" begin
    
    # test collection of OCP in problems/ 
    include("./ci/test_all_ocp.jl") #ok except 2 cases with adnlp_backend=:manual

    # test discretization options
    include("./ci/test_discretization.jl") #ok

    # test NLP modeler / solver options
    include("./ci/test_modeler_solver.jl") #ok except 1 case with adnlp_backend=:manual

end
