using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution
using CTModels: CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTParser: CTParser, @def, prefix!
prefix!(:CTModels) # tell CTParser def macro to use CTModels instead of OptimalControl

using NLPModelsIpopt
using ADNLPModels
using SparseMatrixColorings
using Printf

# load examples library
problem_path = pwd() * "/test/problems"
for problem_file in filter(contains(r".jl$"), readdir(problem_path; join=true))
    include(problem_file)
end

# coloring test function
function coloring_test(ocp; order=NaturalOrder(), grid_size=CTDirect.__grid_size(), disc_method=CTDirect.__disc_method())

    # build DOCP
    time_grid = CTDirect.__time_grid()
    docp = CTDirect.DOCP(ocp; grid_size=grid_size, time_grid=time_grid, disc_method=disc_method)

    # build sparsity pattern
    J = CTDirect.DOCP_Jacobian_pattern(docp)
    H = CTDirect.DOCP_Hessian_pattern(docp)

    ## Coloring for Jacobians
    problem_J = ColoringProblem(; structure=:nonsymmetric, partition=:column)
    order_J = order
    algo_J = GreedyColoringAlgorithm(order_J; decompression=:direct)
    result_J = coloring(J, problem_J, algo_J)
    num_colors_J = ncolors(result_J)

    ## Coloring for Hessians
    problem_H = ColoringProblem(; structure=:symmetric, partition=:column)
    order_H = order
    algo_H = GreedyColoringAlgorithm(order_H; decompression=:substitution, postprocessing=true)
    result_H = coloring(H, problem_H, algo_H)
    num_colors_H = ncolors(result_H)

    return num_colors_J, num_colors_H
end

# batch testing
function batch_coloring_test(; order=NaturalOrder(), target_list=:default, verbose=1, grid_size=CTDirect.__grid_size(), disc_method=CTDirect.__disc_method())

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
        (num_J, num_H) = coloring_test(problem.ocp; order=order, grid_size=grid_size, disc_method=disc_method)
        push!(num_J_list, num_J)
        push!(num_H_list, num_H)
        verbose > 1 && @printf("%-30s J colors %2d    H colors %2d\n", problem.name, num_J, num_H)
    end
    return sum(num_J_list), sum(num_H_list)
end
