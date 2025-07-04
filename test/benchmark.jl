# Benchmark and profiling
using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution
using CTModels: CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTParser: CTParser, @def, prefix!
prefix!(:CTModels) # tell CTParser def macro to use CTModels instead of OptimalControl

using ADNLPModels
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


# solve list of problems, for a given grid size and other options
# verbose <= 1: no output
# verbose > 1: print summary (iter, obj, time)
# verbose > 2: print NLP iterations also
function bench_problem(problem; verbose=1, nlp_solver, kwargs...)

    if verbose > 2
        display = true
    else
        display = false
    end

    # check (will also precompile)
    sol = solve(problem[:ocp], nlp_solver; init=problem[:init], display=display, kwargs...)
    if !isnothing(problem[:obj]) && !isapprox(objective(sol), problem[:obj], rtol=5e-2)
        verbose > 1 && println("\nObjective mismatch for ", problem[:name], ": ", objective(sol), " instead of ", problem[:obj])
        success = false
    else
        verbose > 1 && @printf("\n%-20s: %4d iter %5.2f obj ", problem[:name], iterations(sol), objective(sol))
        success = true
    end
    iter = iterations(sol)

    # time
    time = @belapsed solve($problem[:ocp], $nlp_solver; init=$problem[:init], display=false, $kwargs...)
    verbose > 1 && @printf("%7.2f s\n", time)
    
    return time, iter, success
end


# perform benchmark
function bench(;verbose=1, 
    target_list=:default,
    grid_size_list=[250, 500, 1000, 2500, 5000],
    nlp_solver=:ipopt, 
    kwargs...)

    # load problems for benchmark
    # Note that problems may vary significantly in convergence times...  
    if target_list == :default
        target_list = ["beam", "double_integrator_mintf", "double_integrator_minenergy", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
    elseif target_list == :lagrange_easy
        target_list = [
        "beam",  
        "double_integrator_minenergy", 
        "fuller", 
        "simple_integrator", 
        "vanderpol"]
    elseif target_list == :lagrange_hard
        target_list = [ 
        "bioreactor_1day", 
        "bioreactor_Ndays", 
        "bolza_freetf",  
        "insurance", 
        "parametric", 
        "robbins"]
    elseif target_list == :all
        target_list = ["algal_bacterial", "beam", "bioreactor_1day", "bioreactor_Ndays", "bolza_freetf", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "insurance", "jackson", "parametric", "robbins", "simple_integrator", "swimmer", "vanderpol"]
    elseif target_list == :hard
        target_list = ["algal_bacterial", "bioreactor_1day", "bioreactor_Ndays", "bolza_freetf", "insurance", "swimmer"]
    end
    verbose > 2 && println("\nProblem list: ", target_list)
    problem_list = []
    for problem_name in target_list
        ocp_data = getfield(Main, Symbol(problem_name))()
        push!(problem_list, ocp_data)
    end

    # solve problem list for all grid sizes
    verbose > 2 && println("\nGrid size list: ", grid_size_list)
    t_bench = []
    i_bench = []
    s_bench = []
    for problem in problem_list
        verbose > 1 && @printf("\nTesting problem %-20s for grid size ", problem[:name])
        t_list = []
        i_list = []
        s_list = []
        for grid_size in grid_size_list
            verbose > 1 && @printf("%d ", grid_size)
            time, iter, success = bench_problem(problem; grid_size=grid_size, verbose=verbose-1, nlp_solver=nlp_solver, kwargs...)
            append!(t_list, time)
            append!(i_list, iter)
            append!(s_list, success)
        end
        push!(t_bench, t_list)
        push!(i_bench, i_list)
        push!(s_bench, s_list)
    end

    # display: 1 row per problem with (t,i,s) for each grid size in columns
    # plus last row for total over problem set
    if verbose > 0
        i = 1
        for problem in problem_list
            @printf("\n%-20s ", problem[:name])
            for j=1:length(grid_size_list)
                if s_bench[i][j]
                    @printf("%4.1f (%4d) ", t_bench[i][j], i_bench[i][j])
                else
                    @printf("FAIL     ")
                end
            end
            i = i + 1
        end
    end
    #+++total, quid des fail ?
    println(t_bench)
    println(i_bench)
    println(i_bench[:][1])
    @printf("\nBENCHMARK             ")
    #for j=1:length(grid_size_list)
    #    @printf("%4.1f (%4d)   ", sum(t_bench[:][j]), sum(i_bench[:][j]))
    #end
    return
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
