# Benchmark and profiling
using CTBase #? still needed ?
using CTParser: CTParser, @def
using CTModels:
    CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution

using ADNLPModels
using NLPModelsIpopt
using MadNLPMumps

using LinearAlgebra
using Printf
using Plots

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

# check a specific example
function check_problem(prob; kwargs...)
    sol = solve(prob.ocp; init=prob.init, kwargs...)
    @test sol.objective ≈ prob.obj rtol = 1e-2
end

# tests to check allocations in particular
function init(ocp; grid_size, disc_method)
    docp = CTDirect.DOCP(
        ocp; grid_size=grid_size, time_grid=CTDirect.__time_grid(), disc_method=disc_method
    )
    xu = CTDirect.DOCP_initial_guess(docp)
    return docp, xu
end

function test_unit(
    ocp;
    test_obj=true,
    test_cons=true,
    test_trans=true,
    test_solve=true,
    warntype=false,
    jet=false,
    profile=false,
    grid_size=100,
    disc_method=:trapeze,
)
    if profile
        Profile.Allocs.clear()
    end

    # define problem and variables
    docp, xu = init(ocp; grid_size=grid_size, disc_method=disc_method)
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
        @btime direct_transcription(
            $ocp, grid_size=($grid_size), disc_method=($disc_method)
        )
    end

    # solve
    if test_solve
        print("Solve");
        @btime solve(
            $ocp, display=false, grid_size=($grid_size), disc_method=($disc_method)
        )
    end

    return nothing
end

# solve given problem, return convergence data and solution
# verbose <= 1: no output
# verbose > 1: print summary (iter, obj, time)
# verbose > 2: print NLP iterations also
function bench_problem(problem; verbose=1, nlp_solver, grid_size, kwargs...)
    if verbose > 2
        display = true
    else
        display = false
    end

    # check (will also precompile)
    time = @elapsed sol = solve(
        problem[:ocp],
        nlp_solver;
        init=problem[:init],
        display=display,
        grid_size=grid_size,
        kwargs...,
    )
    if !CTModels.successful(sol) ||
        (!isnothing(problem[:obj]) && !isapprox(objective(sol), problem[:obj]; rtol=5e-2))
        success = false
        iter = min(iterations(sol), 999) # to fit 3-digit print 
        println(
            "\nFailed ",
            problem[:name],
            " for grid size ",
            grid_size,
            " at iter ",
            iter,
            " obj ",
            objective(sol),
            " vs ",
            problem[:obj],
        )
    else
        success = true
        iter = iterations(sol)
        verbose > 1 && @printf(
            "%-20s: %4d iter %5.2f obj ",
            problem[:name],
            iterations(sol),
            objective(sol)
        )
        # time
        time = @belapsed solve(
            $problem[:ocp],
            $nlp_solver;
            init=$problem[:init],
            display=false,
            grid_size=($grid_size),
            $kwargs...,
        )
        verbose > 1 && @printf("%7.2f s\n", time)
    end

    return time, iter, success, sol
end

# perform benchmark
function bench(;
    verbose=1,
    target_list=:all,
    grid_size_list=[250, 500, 1000],
    nlp_solver=:ipopt,
    return_sols=false,
    save_sols=false,
    kwargs...,
)

    # load problems for benchmark
    # Note that problems may vary significantly in convergence times...  
    if target_list == :lagrange_easy
        target_list = [
            "beam",
            "double_integrator_minenergy",
            "fuller",
            "simple_integrator",
            "vanderpol",
        ]
    elseif target_list == :lagrange_hard
        target_list = [
            "bioreactor_1day",
            "bioreactor_Ndays",
            "bolza_freetf",
            "insurance", #only converge when final control is present (mixed path constraint) 
            "parametric",
            "robbins",
        ]
    elseif target_list == :lagrange_all
        target_list = [
            "beam",
            "bioreactor_1day",
            "bioreactor_Ndays",
            "bolza_freetf",
            "double_integrator_e",
            "fuller",
            "parametric",
            "robbins",
            "simple_integrator",
            "vanderpol",
        ]
    elseif target_list == :hard
        target_list = [
            "action",
            "glider",
            "moonlander",
            "quadrotor",
            "schlogl",
            "space_shuttle",
            "truck_trailer",
        ]

    elseif target_list == :all
        target_list = [
            "algal_bacterial",
            "beam",
            "bioreactor_1day",
            "bioreactor_Ndays",
            "bolza_freetf",
            "double_integrator_mintf",
            "double_integrator_minenergy",
            "double_integrator_freet0tf",
            "fuller",
            "goddard",
            "goddard_all",
            #"insurance", fail unless final control
            "jackson",
            "parametric",
            "robbins",
            "simple_integrator",
            #"swimmer", #much slower then others #fail for madnlpmumps
            "vanderpol",
        ]
    elseif target_list == :hard
        target_list = [
            "algal_bacterial",
            "bioreactor_1day",
            "bioreactor_Ndays",
            "bolza_freetf",
            "insurance",
            "swimmer",
        ]
    end
    verbose > 1 && println("Problem list: ", target_list)
    problem_list = []
    for problem_name in target_list
        ocp_data = getfield(Main, Symbol(problem_name))()
        push!(problem_list, ocp_data)
    end

    # solve problem list for all grid sizes
    verbose > 1 && println("Grid size list: ", grid_size_list)
    t_bench = zeros(Float64, (length(problem_list), length(grid_size_list)))
    i_bench = zeros(Int, (length(problem_list), length(grid_size_list)))
    s_bench = zeros(Bool, (length(problem_list), length(grid_size_list)))
    solutions = Array{Any}(undef, (length(problem_list), length(grid_size_list)))
    i = 1
    for problem in problem_list
        verbose > 1 && @printf("Testing problem %-22s for grid size ", problem[:name])
        j = 1
        for grid_size in grid_size_list
            verbose > 1 && @printf("%d ", grid_size)
            flush(stdout)
            time, iter, success, sol = bench_problem(
                problem;
                grid_size=grid_size,
                verbose=verbose-1,
                nlp_solver=nlp_solver,
                kwargs...,
            )
            t_bench[i, j] = time
            i_bench[i, j] = iter
            s_bench[i, j] = success
            solutions[i, j] = sol
            j = j + 1
        end
        verbose > 1 && println("")
        i = i + 1
    end

    # display: 1 row per problem with (t,i,s) for each grid size in columns
    # plus last row for total over problem set
    if verbose > 0
        i = 1
        for problem in problem_list
            @printf("%-22s", problem[:name])
            for j in 1:length(grid_size_list)
                if s_bench[i, j]
                    @printf("%6.2f(%3d) ", t_bench[i, j], i_bench[i, j])
                else
                    @printf("  FAIL(%3d) ", i_bench[i, j])
                end
            end
            println("")
            i = i + 1
        end
    end

    # summary
    @printf("SUCCESS %2d/%2d    ", sum(s_bench), length(s_bench))
    for j in 1:length(grid_size_list)
        @printf("%6.2f(%3d) ", sum(t_bench[:, j]), sum(i_bench[:, j]))
    end
    println("")

    if return_sols
        return solutions
    else
        return nothing
    end
end

# custom bench calls
function bench_custom()
    disc_list = [
        #:euler,
        #:euler_implicit,
        :trapeze,
        :midpoint,
        #:gauss_legendre_2,
        #:gauss_legendre_3
    ]

    target_list = ["goddard"] #:hard
    grid_size_list=[250] #, 500, 1000]
    verbose = 1

    solutions = Dict{Symbol,Any}()

    for disc in disc_list
        lagrange_to_mayer=true
        @printf("Bench %s / %s Lag2Mayer ", target_list, disc)
        println(lagrange_to_mayer, " Grid ", grid_size_list)
        solutions[disc] = bench(;
            target_list=target_list,
            grid_size_list=grid_size_list,
            disc_method=disc,
            verbose=verbose,
            lagrange_to_mayer=lagrange_to_mayer,
        )
        flush(stdout)
        println("")

        # plot
        if disc == disc_list[1]
            p = plot(solutions[disc][1, 1], :control; label=String(disc))
        else
            p = plot!(solutions[disc][1, 1], :control; label=String(disc))
        end
        display(p)

        #=
        lagrange_to_mayer=false
        @printf("Bench %s / %s Lag2Mayer ", target_list, disc)
        println(lagrange_to_mayer, " Grid ", grid_size_list)
        sol2 = bench(target_list=target_list, grid_size_list=grid_size_list, disc_method=disc, verbose=verbose, lagrange_to_mayer=lagrange_to_mayer)
        flush(stdout)
        println("")
        =#
    end
end
