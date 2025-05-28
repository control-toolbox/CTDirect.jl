using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution
using CTModels: CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTParser: CTParser, @def, prefix!
prefix!(:CTModels) # tell CTParser def macro to use CTModels instead of OptimalControl

using NLPModelsIpopt
using ADNLPModels
using SparseMatrixColorings
using Printf

using CliqueTrees

# load examples library
problem_path = pwd() * "/test/problems"
for problem_file in filter(contains(r".jl$"), readdir(problem_path; join = true))
    include(problem_file)
end

function get_patterns(ocp; grid_size=CTDirect.__grid_size(), disc_method=CTDirect.__disc_method(), adnlp_backend=CTDirect.__adnlp_backend())

    docp, nlp = direct_transcription(ocp; grid_size=grid_size, disc_method=disc_method, adnlp_backend=adnlp_backend)
    J = get_sparsity_pattern(nlp, :jacobian)
    H = get_sparsity_pattern(nlp, :hessian) + transpose(get_sparsity_pattern(nlp, :hessian))
   
    return J, H 
end

# coloring test functions
function test_Jacobian_coloring(J, order)

    problem_J = ColoringProblem(; structure=:nonsymmetric, partition=:column)
    if order == PerfectEliminationOrder()
        order_J = NaturalOrder() 
    else 
        order_J = order
    end
    algo_J = GreedyColoringAlgorithm(order_J; decompression=:direct)
   
    return ncolors(coloring(J, problem_J, algo_J))
end

function test_Hessian_coloring(H, order)
   
    problem_H = ColoringProblem(; structure=:symmetric, partition=:column)
    order_H = order
    algo_H = GreedyColoringAlgorithm(order_H; decompression=:substitution, postprocessing=true) # check vs direct ?
   
    return ncolors(coloring(H, problem_H, algo_H))
end

# batch testing
# possible orders: NaturalOrder(), RandomOrder(), LargestFirst(), SmallestLast(), IncidenceDegree(), DynamicLargestFirst(), PerfectEliminationOrder() (for Hessian with substitution decomposition)
function batch_coloring_test(; order = NaturalOrder(), target_list = :default, verbose = 1, grid_size = CTDirect.__grid_size(), disc_method = CTDirect.__disc_method(), adnlp_backend = CTDirect.__adnlp_backend())

    if target_list == :default
        target_list = ["beam", "double_integrator_mintf", "double_integrator_minenergy", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
    end

    verbose > 1 && println("\nProblem list: ", target_list)
    problem_list = []
    for problem_name in target_list
        ocp_data = getfield(Main, Symbol(problem_name))()
        push!(problem_list, ocp_data)
    end

    num_J_list = []
    num_H_list = []
    for problem in problem_list
        J, H = get_patterns(problem.ocp; grid_size=grid_size, disc_method=disc_method, adnlp_backend=adnlp_backend)
        num_J = test_Jacobian_coloring(J, order)
        num_H = test_Hessian_coloring(H, order)
        push!(num_J_list, num_J)
        push!(num_H_list, num_H)
        verbose > 1 && @printf("%-30s J colors %2d    H colors %2d\n", problem.name, num_J, num_H)
    end
    return sum(num_J_list), sum(num_H_list)
end

