# Benchmark
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt

using Printf

using MKL # Replace OpenBLAS with Intel MKL +++ should be an option

#using MadNLPMumps

#######################################################
# load examples library
problem_path = pwd() * "/test/problems"
for problem_file in filter(contains(r".jl$"), readdir(problem_path; join = true))
    include(problem_file)
end

function bench(;
    nlp_solver = :ipopt,
    linear_solver = nothing,
    tol = 1e-8,
    grid_size = 1000,
    precompile = true,
    display = false,
    verbose = true,
    disc_method = :trapeze
)

    #######################################################
    # set (non) linear solvers and backends

    # linear solver for ipopt: default mumps; spral, ma27, ma57, ma77, ma86, ma97
    if nlp_solver == :ipopt && isnothing(linear_solver)
        linear_solver = "mumps"
    end

    # linear solver for madnlp: default umfpack; MumpsSolver
    if nlp_solver == :madnlp && isnothing(linear_solver)
        linear_solver = "UmfpackSolver"
    end

    verbose && @printf("Profile: NLP Solver %s with linear solver %s\n", nlp_solver, linear_solver)

    # blas backend (cf using MKL above, should be option...)
    verbose && @printf("Blas config: %s\n", LinearAlgebra.BLAS.lbt_get_config())

    # settings
    verbose &&
        @printf("Settings: tol=%g grid_size=%d precompile=%s\n\n", tol, grid_size, precompile)

    # load problems for benchmark
    names_list = ["beam", "bioreactor_1day", "fuller", "goddard", "jackson", "vanderpol"]
    problem_list = []
    for problem_name in names_list
        ocp_data = getfield(Main, Symbol(problem_name))()
        push!(problem_list, ocp_data)
    end

    #######################################################
    # precompile if required
    if precompile
        t_precomp = 0.0
        verbose && print("Precompilation: ")
        for problem in problem_list
            verbose && @printf("%s ", problem[:name])
            t = @elapsed direct_solve(
                problem[:ocp],
                nlp_solver,
                linear_solver = linear_solver,
                max_iter = 0,
                display = display,
                disc_method = disc_method
            )
            t_precomp += t
        end
        verbose && @printf("\nPrecompilation total time %6.2f\n\n", t_precomp)
    end

    #######################################################
    # solve examples with timer and objective check
    t_list = []
    for problem in problem_list
        t = @elapsed sol = direct_solve(
            problem[:ocp],
            nlp_solver,
            init = problem[:init],
            display = display,
            linear_solver = linear_solver,
            grid_size = grid_size,
            tol = tol,
            disc_method = disc_method
        )
        if !isnothing(problem[:obj]) && !isapprox(sol.objective, problem[:obj], rtol = 5e-2)
            error(
                "Objective mismatch for ",
                problem[:name],
                ": ",
                sol.objective,
                " instead of ",
                problem[:obj],
            )
        else
            verbose && @printf(
                "%-30s completed in %6.2f s after %4d iterations\n",
                problem[:name],
                t,
                sol.iterations
            )
            append!(t_list, t)
        end
    end

    #######################################################
    # print total time
    total_time = sum(t_list)
    verbose && @printf("\nTotal time (s): %6.2f\n", total_time)

    # return also full text ouptut ?
    return total_time
end

# +++ put repeat directly in bench()
function bench_average(; repeat = 4, verbose = false, kwargs...)

    # execute series of benchmark runs
    t_list = []
    for i = 1:repeat
        t = bench(; verbose = verbose, kwargs...)
        append!(t_list, t)
        verbose && @printf("Run %d / %d: time (s) = %6.2f\n", i, repeat, t)
    end

    # print / return average total time
    avg_time = sum(t_list) / length(t_list)
    verbose && @printf("Average time (s): %6.2f\n", avg_time)
    return avg_time
end

function bench_series(; grid_size_list = [250, 500, 1000, 2500, 5000], kwargs...)
    println(grid_size_list)
    t_list = []
    for grid_size in grid_size_list
        t = bench_average(; grid_size = grid_size, kwargs...)
        append!(t_list, t)
        @printf("Grid size %d: time (s) = %6.1f\n", grid_size, t)
    end
    #return t_list
end
