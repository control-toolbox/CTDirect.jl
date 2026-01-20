# solve tests

include("../problems/beam.jl")

function test_ctdirect_solve()

    # common options
    ipopt_options = Dict(
        :max_iter => 1000,
        :tol => 1e-6,
        :print_level => 3,
        :mu_strategy => "adaptive",
        :linear_solver => "Mumps",
        :sb => "yes",
    )

    # ========================================================================
    # INTEGRATION: Direct beam OCP with Collocation (Ipopt pieces)
    # ========================================================================
    Test.@testset "integration: beam_docp" verbose=VERBOSE showtiming=SHOWTIMING begin

        # load and discretize OCP
        beam_data = beam()
        ocp = beam_data.ocp
        init = CTModels.initial_guess(ocp; beam_data.init...)
        discretizer = CTDirect.Collocation()
        docp = CTDirect.discretize(ocp, discretizer)
        Test.@test docp isa CTModels.DiscretizedOptimalControlProblem

        # NLP modelers
        modelers = [CTModels.ADNLPModeler()]#, CTModels.ExaModeler()]
        modelers_names = ["ADNLPModeler"] #, "ExaModeler (CPU)"]

        # solve NLP with ipopt specific solver and NLP modelers
        Test.@testset "NLP level (solve_with_ipopt)" verbose=VERBOSE showtiming=SHOWTIMING begin
            for (modeler, modeler_name) in zip(modelers, modelers_names)
                Test.@testset "$(modeler_name)" verbose=VERBOSE showtiming=SHOWTIMING begin
                    nlp = CTModels.nlp_model(docp, init, modeler)
                    stats = CTSolvers.solve_with_ipopt(nlp; ipopt_options...)
                    sol = CTModels.ocp_solution(docp, stats, modeler)
                    Test.@test sol isa CTModels.Solution
                    Test.@test CTModels.successful(sol)
                    Test.@test isfinite(CTModels.objective(sol))
                    Test.@test CTModels.objective(sol) ≈ beam_data.obj atol=1e-2
                end
            end
        end

        # solve DOCP with common solve and NLP modelers
        Test.@testset "DOCP level (solve)" verbose=VERBOSE showtiming=SHOWTIMING begin
            for (modeler, modeler_name) in zip(modelers, modelers_names)
                Test.@testset "$(modeler_name)" verbose=VERBOSE showtiming=SHOWTIMING begin
                    solver = CTSolvers.IpoptSolver(; ipopt_options...)
                    sol = CommonSolve.solve(docp, init, modeler, solver; display=false)
                    Test.@test sol isa CTModels.Solution
                    Test.@test CTModels.successful(sol)
                    Test.@test isfinite(CTModels.objective(sol))
                    Test.@test CTModels.objective(sol) ≈ beam_data.obj atol=1e-2
                    Test.@test CTModels.iterations(sol) <= ipopt_options[:max_iter]
                    Test.@test CTModels.constraints_violation(sol) <= 1e-6
                end
            end
        end

        # check_problem toplevel function
        Test.@testset "DOCP level (check_problem)" verbose=VERBOSE showtiming=SHOWTIMING begin
            check_problem(beam(); display=true)
        end
    end



    # ========================================================================
    # INTEGRATION: Direct Goddard OCP with Collocation (Ipopt pieces)
    #= ========================================================================
    Test.@testset "integration: goddard_docp" verbose=VERBOSE showtiming=SHOWTIMING begin
        if !isdefined(Main, :goddard)
            include("./problems/goddard.jl")
        end
        gdata = goddard()
        ocp_g = gdata.ocp
        init_g = CTSolvers.initial_guess(ocp_g; gdata.init...)
        discretizer_g = CTSolvers.Collocation()
        docp_g = CTSolvers.discretize(ocp_g, discretizer_g)

        Test.@test docp_g isa CTSolvers.DiscretizedOptimalControlProblem

        modelers = [CTSolvers.ADNLPModeler(), CTSolvers.ExaModeler()]
        modelers_names = ["ADNLPModeler", "ExaModeler (CPU)"]

        # ocp_solution from DOCP using solve_with_ipopt
        Test.@testset "ocp_solution from DOCP (Ipopt)" verbose=VERBOSE showtiming=SHOWTIMING begin
            for (modeler, modeler_name) in zip(modelers, modelers_names)
                Test.@testset "$(modeler_name)" verbose=VERBOSE showtiming=SHOWTIMING begin
                    nlp = CTModels.nlp_model(docp_g, init_g, modeler)
                    stats = CTSolvers.solve_with_ipopt(nlp; ipopt_options...)
                    sol = CTSolvers.ocp_solution(docp_g, stats, modeler)
                    Test.@test sol isa CTModels.Solution
                    Test.@test CTModels.successful(sol)
                    Test.@test isfinite(CTModels.objective(sol))
                    Test.@test CTModels.objective(sol) ≈ gdata.obj atol=1e-4
                end
            end
        end

        # DOCP level: CommonSolve.solve(docp_g, init_g, modeler, solver)
        Test.@testset "DOCP level (Ipopt)" verbose=VERBOSE showtiming=SHOWTIMING begin
            for (modeler, modeler_name) in zip(modelers, modelers_names)
                Test.@testset "$(modeler_name)" verbose=VERBOSE showtiming=SHOWTIMING begin
                    solver = CTSolvers.IpoptSolver(; ipopt_options...)
                    sol = CommonSolve.solve(docp_g, init_g, modeler, solver; display=false)
                    Test.@test sol isa CTModels.Solution
                    Test.@test CTModels.successful(sol)
                    Test.@test isfinite(CTModels.objective(sol))
                    Test.@test CTModels.objective(sol) ≈ gdata.obj atol=1e-4
                    Test.@test CTModels.iterations(sol) <= ipopt_options[:max_iter]
                    Test.@test CTModels.constraints_violation(sol) <= 1e-6
                end
            end
        end
    end=#
end