# Benchmark
using CTDirect
using CTBase

using LinearAlgebra
using NLPModelsIpopt

using MKL # Replace OpenBLAS with Intel MKL +++ should be an option

using BenchmarkTools
using Printf


#using MadNLPMumps

#######################################################
# load examples library
problem_path = pwd() * "/test/problems"
for problem_file in filter(contains(r".jl$"), readdir(problem_path; join = true))
    include(problem_file)
end



function bench_list(problem_list; verbose=2, nlp_solver, linear_solver, kwargs...)

    #######################################################
    # solve examples with timer and objective check
    t_list = []
    for problem in problem_list

        # check
        sol = direct_solve(problem[:ocp], nlp_solver; init=problem[:init], display=false, kwargs...)
        if !isnothing(problem[:obj]) && !isapprox(sol.objective, problem[:obj], rtol = 5e-2)
            error("Objective mismatch for ",problem[:name],": ",sol.objective," instead of ",problem[:obj])
        else
            verbose > 1 && @printf("%-30s: %4d iter ", problem[:name], sol.iterations)
        end

        # time
        t = @belapsed direct_solve($problem[:ocp], $nlp_solver; init=$problem[:init], display=false, $kwargs...)
        append!(t_list, t)
        verbose > 1 && @printf("%7.2f s\n", t)
    end

    return sum(t_list)
end


function bench(;grid_size_list = [250, 500, 1000, 2500, 5000], verbose = 1, nlp_solver=:ipopt, linear_solver=nothing, names_list = :default, kwargs...)

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
    verbose > 1 && @printf("Profile: NLP Solver %s with linear solver %s\n", nlp_solver, linear_solver)
    # blas backend (cf using MKL above, should be option...)
    verbose > 1 && @printf("Blas config: %s\n", LinearAlgebra.BLAS.lbt_get_config())

    # load problems for benchmark
    if names_list == :default
        names_list = ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "robbins", "simple_integrator", "vanderpol"]
    elseif names_list == :quick
        names_list = ["beam", "double_integrator_mintf", "fuller", "jackson", "robbins", "simple_integrator", "vanderpol"]
    elseif names_list == :all 
        names_list = ["algal_bacterial", "beam", "bioreactor_1day", "bioreactor_Ndays", "bolza_freetf", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "insurance", "jackson", "robbins", "simple_integrator", "swimmer", "vanderpol"]
    elseif names_list == :hard
        names_list = ["algal_bacterial", "bioreactor_1day", "bioreactor_Ndays", "bolza_freetf", "goddard_all", "insurance", "swimmer"]
    end
    println("Problem list: ", names_list)
    problem_list = []
    for problem_name in names_list
        ocp_data = getfield(Main, Symbol(problem_name))()
        push!(problem_list, ocp_data)
    end

    println("Grid size list: ", grid_size_list)
    t_list = []
    for grid_size in grid_size_list
        t = bench_list(problem_list; grid_size=grid_size, verbose=verbose, nlp_solver=nlp_solver, linear_solver=linear_solver, kwargs...)
        append!(t_list, t)
        @printf("Grid size %d: time (s) = %6.1f\n", grid_size, t)
    end

end
