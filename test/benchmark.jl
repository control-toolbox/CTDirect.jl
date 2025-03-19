# Benchmark
using CTDirect
import CTModels

using LinearAlgebra
using MadNLP
using NLPModelsIpopt

using MKL # Replace OpenBLAS with Intel MKL +++ should be an option

using BenchmarkTools
#using Plots
using Printf
using Profile
using PProf
using JET
using Test


#######################################################
# load examples library
problem_path = pwd() * "/test/problems"
for problem_file in filter(contains(r".jl$"), readdir(problem_path; join = true))
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
        sol = direct_solve(problem[:ocp], nlp_solver; init=problem[:init], display=display, kwargs...)
        if !isnothing(problem[:obj]) && !isapprox(sol.objective, problem[:obj], rtol = 5e-2)
            error("Objective mismatch for ",problem[:name],": ",sol.objective," instead of ",problem[:obj])
        else
            verbose > 2 && @printf("%-30s: %4d iter ", problem[:name], sol.iterations)
        end

        # time
        t = @belapsed direct_solve($problem[:ocp], $nlp_solver; init=$problem[:init], display=false, $kwargs...)
        append!(t_list, t)
        verbose > 2 && @printf("%7.2f s\n", t)
    end

    return sum(t_list)
end


function bench(; grid_size_list = [250, 500, 1000, 2500, 5000], verbose = 1, nlp_solver=:ipopt, linear_solver=nothing, names_list = :default, kwargs...)

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
    if names_list == :default
        names_list = ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
    elseif names_list == :quick
        names_list = ["beam", "double_integrator_mintf", "fuller", "jackson", "simple_integrator", "vanderpol"]
    elseif names_list == :all 
        names_list = ["algal_bacterial", "beam", "bioreactor_1day", "bioreactor_Ndays", "bolza_freetf", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "insurance", "jackson", "robbins", "simple_integrator", "swimmer", "vanderpol"]
    elseif names_list == :hard
        names_list = ["algal_bacterial", "bioreactor_1day", "bioreactor_Ndays", "bolza_freetf", "goddard_all", "insurance", "swimmer"]
    end
    verbose > 1 && println("Problem list: ", names_list)
    problem_list = []
    for problem_name in names_list
        ocp_data = getfield(Main, Symbol(problem_name))()
        push!(problem_list, ocp_data)
    end

    verbose > 1 && println("Grid size list: ", grid_size_list)
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
    #= OK, same as calling the functions with docp
    NLP_objective = (xu) -> CTDirect.DOCP_objective(xu, docp)
    NLP_constraints! = (c, xu) -> CTDirect.DOCP_constraints!(c, xu, docp) =#
    c = fill(666.666, docp.dim_NLP_constraints)
    work = similar(xu, docp.dim_NLP_x)

    # DOCP_objective
    if test_obj
        print("Objective"); @btime CTDirect.DOCP_objective($xu, $docp)
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
        print("Constraints"); @btime CTDirect.DOCP_constraints!($c, $xu, $docp)
        any(c.==666.666) && error("undefined values in constraints ",c)
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
        print("Transcription"); @btime direct_transcription($ocp, grid_size=$grid_size, disc_method=$disc_method)
    end

    # solve
    if test_solve
        sol = direct_solve(ocp, display=false, grid_size=grid_size, disc_method=disc_method)
        if !isapprox(sol.objective, prob.obj, rtol=1e-2)
            error("objective mismatch: ", sol.objective, " vs ", prob.obj)
        end
        print("Solve"); @btime direct_solve($ocp, display=false, grid_size=$grid_size, disc_method=$disc_method)
    end

end

