# Benchmark
include("../test/deps.jl")
using Printf

using MKL # Replace OpenBLAS with Intel MKL +++ should be an option

using MadNLPMumps

function bench(;nlp_solver = :ipopt, linear_solver = nothing,
    tol=1e-8, grid_size=1000, precompile = true, display=false)

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

    @printf("Profile: NLP Solver %s with linear solver %s\n", nlp_solver, linear_solver)

    # blas backend (cf using MKL above, should be option...)
    @printf("Blas config: %s\n", LinearAlgebra.BLAS.lbt_get_config())

    # settings
    @printf("Settings: tol=%g grid_size=%d precompile=%s\n\n", tol, grid_size, precompile)

    #######################################################
    # load examples library
    problem_path = pwd()*"/test/problems"
    for problem_file in filter(contains(r".jl$"), readdir(problem_path; join=true))
        include(problem_file)
    end

    # load problems for benchmark
    names_list = ["beam", "bioreactor_1day", "fuller", "goddard", "jackson", "vanderpol"]
    problem_list = []
    for problem_name in names_list
        ocp_data = getfield(Main, Symbol(problem_name))()
        #@printf("%s ", ocp_data.name)
        push!(problem_list,ocp_data)
    end
    #println("")

    #######################################################
    # precompile if required
    if precompile
        t_precomp = 0.
        print("Precompilation: ")
        for problem in problem_list
            @printf("%s ",problem[:name])
            t = @elapsed solve(problem[:ocp], nlp_solver, linear_solver=linear_solver, max_iter=0, display=display)
            t_precomp += t
        end
        @printf("\nPrecompilation total time %6.2f\n\n",t_precomp)
    end

    #######################################################
    # solve examples with timer and objective check
    t_list = []
    #println("Benchmark:")
    for problem in problem_list
        t = @elapsed sol = solve(problem[:ocp], nlp_solver, init=problem[:init], display=display, linear_solver=linear_solver, grid_size=grid_size, tol=tol)
        if !isnothing(problem[:obj]) && !isapprox(sol.objective, problem[:obj], rtol=5e-2)
            error("Objective mismatch for ", problem[:name], ": ", sol.objective, " instead of ", problem[:obj])
        else
            @printf("%-30s completed in %6.2f s after %4d iterations\n",problem[:name],t,sol.iterations)
            append!(t_list,t)
        end
    end

    #######################################################
    # print total time
    total_time = sum(t_list)
    @printf("\nTotal time (s): %6.2f", total_time)

    # return also full text ouptut ?
    return total_time

end
