# Benchmark and profiling
using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution
using CTModels: CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTParser: CTParser, @def, set_prefix
set_prefix(:CTModels) # tell CTParser def macro to use CTModels instead of OptimalControl

using NLPModelsIpopt

using LinearAlgebra
using Printf

using BenchmarkTools
using JET
using Profile
using PProf
using Test # to run individual test scripts if needed


#######################################################
# load examples library
problem_path = pwd() * "/test/problems"
for problem_file in filter(contains(r".jl$"), readdir(problem_path; join=true))
    include(problem_file)
end

function bench_list(problem_list; verbose=1, nlp_solver, linear_solver, kwargs...)

    if verbose > 3
        display = true
    else
        display = false
    end

    # solve examples with timer and objective check
    t_list = []
    for problem in problem_list

        # check (will also precompile)
        sol = solve(problem[:ocp], nlp_solver; init=problem[:init], display=display, kwargs...)
        if !isnothing(problem[:obj]) && !isapprox(objective(sol), problem[:obj], rtol=5e-2)
            error("Objective mismatch for ", problem[:name], ": ", objective(sol), " instead of ", problem[:obj])
        else
            verbose > 2 && @printf("%-30s: %4d iter %5.2f obj ", problem[:name], iterations(sol), objective(sol))
        end

        # time
        t = @belapsed solve($problem[:ocp], $nlp_solver; init=$problem[:init], display=false, $kwargs...)
        append!(t_list, t)
        verbose > 2 && @printf("%7.2f s\n", t)
    end

    return sum(t_list)
end


function bench(; grid_size_list=[250, 500, 1000, 2500, 5000], verbose=1, nlp_solver=:ipopt, linear_solver=nothing, target_list=:default, kwargs...)

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
    # Note that problems may vary significantly in convergence times...  
    if target_list == :default
        target_list = ["beam", "double_integrator_mintf", "double_integrator_minenergy", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
    elseif target_list == :quick
        target_list = ["beam", "double_integrator_mintf", "fuller", "jackson", "simple_integrator", "vanderpol"]
    elseif target_list == :all
        target_list = ["algal_bacterial", "beam", "bioreactor_1day", "bioreactor_Ndays", "bolza_freetf", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "insurance", "jackson", "parametric", "robbins", "simple_integrator", "swimmer", "vanderpol"]
    elseif target_list == :hard
        target_list = ["algal_bacterial", "bioreactor_1day", "bioreactor_Ndays", "bolza_freetf", "insurance", "swimmer"]
    end
    verbose > 1 && println("\nProblem list: ", target_list)
    problem_list = []
    for problem_name in target_list
        ocp_data = getfield(Main, Symbol(problem_name))()
        push!(problem_list, ocp_data)
    end

    verbose > 1 && println("\nGrid size list: ", grid_size_list)
    t_list = []
    for grid_size in grid_size_list
        t = bench_list(problem_list; grid_size=grid_size, verbose=verbose, nlp_solver=nlp_solver, linear_solver=linear_solver, kwargs...)
        append!(t_list, t)
        @printf("Grid size %6d: time (s) = %6.1f\n", grid_size, t)
    end
    #return t_list
end


# tests to check allocations in particular
function init(ocp; grid_size, disc_method)
    docp = CTDirect.DOCP(ocp, grid_size=grid_size, time_grid=CTDirect.__time_grid(), disc_method=disc_method)
    xu = CTDirect.DOCP_initial_guess(docp)
    return docp, xu
end


function test_unit(ocp; test_obj=true, test_cons=true, test_trans=true, test_solve=true, warntype=false, jet=false, profile=false, grid_size=100, disc_method=:trapeze)

    if profile
        Profile.Allocs.clear()
    end

    # define problem and variables
    docp, xu = init(ocp, grid_size=grid_size, disc_method=disc_method)
    disc = docp.discretization
    c = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dims.NLP_x)

    # DOCP_objective
    if test_obj
        print("Objective");
        @btime CTDirect.DOCP_objective($xu, $docp)
        warntype && @code_warntype CTDirect.DOCP_objective(xu, docp)
        jet && display(@report_opt CTDirect.DOCP_objective(xu, docp))
        if profile
            Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_objective(xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    # DOCP_constraints
    if test_cons
        print("Constraints");
        @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        any(c .== 666.666) && error("undefined values in constraints ", c)
        warntype && @code_warntype CTDirect.DOCP_constraints!(c, xu, docp)
        jet && display(@report_opt CTDirect.DOCP_constraints!(c, xu, docp))
        if profile
            Profile.Allocs.@profile sample_rate=1.0 CTDirect.DOCP_constraints!(c, xu, docp)
            results = Profile.Allocs.fetch()
            PProf.Allocs.pprof()
        end
    end

    # transcription
    if test_trans
        print("Transcription");
        @btime direct_transcription($ocp, grid_size=($grid_size), disc_method=($disc_method))
    end

    # solve
    if test_solve
        print("Solve");
        @btime solve($ocp, display=false, grid_size=($grid_size), disc_method=($disc_method))
    end

    return nothing
end
