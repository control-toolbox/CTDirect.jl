# runtests.jl
using Test
include("test_common.jl")

function check_problem(prob; 
    modeler = CTModels.ADNLPModeler(),
    solver = CTSolvers.IpoptSolver(; ipopt_options...),
    display=false)

    sol = solve_problem(prob; modeler, solver, display)
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
    #include("./ci/test_core_types.jl")
    #test_ctdirect_core_types()
    #include("./ci/test_discretization_api.jl")
    #test_ctdirect_discretization_api()
    #include("./ci/test_collocation.jl")
    #test_ctdirect_collocation()

    include("./ci/test_solve.jl")
    test_solve()
    
    #include("./ci/test_all_ocp.jl")
end
