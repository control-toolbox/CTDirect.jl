# solve tests

include("../problems/goddard.jl")

function test_solve()

    # common options
    max_iter = 1000
    tol = 1e-6
    ipopt_options = Dict(
        :max_iter => max_iter,
        :tol => tol,
        :print_level => 3,
        :mu_strategy => "adaptive",
        :linear_solver => "Mumps",
        :sb => "yes",
    )

    # ========================================================================
    # INTEGRATION: Direct beam OCP with Collocation (Ipopt pieces)
    # ========================================================================
    Test.@testset "integration: beam_docp" verbose=VERBOSE showtiming=SHOWTIMING begin

        #= load and discretize OCP
        beam_data = beam2() # use exa-compatible version
        ocp = beam_data.ocp
        init = CTModels.initial_guess(ocp; beam_data.init...)
        discretizer = CTDirect.Collocation()
        docp = CTDirect.discretize(ocp, discretizer)
        Test.@test docp isa CTModels.DiscretizedModel

        # NLP solver
        solvers = [CTSolvers.IpoptSolver(; ipopt_options...)]
        solvers_names = ["Ipopt"]

        # NLP modelers
        modelers = [CTModels.ADNLP(), CTModels.ExaModeler()]
        modelers_names = ["ADNLP", "ExaModeler (CPU)"]

        # solve DOCP with common solve and NLP modelers
        Test.@testset "DOCP level (solve)" verbose=VERBOSE showtiming=SHOWTIMING begin
            for (modeler, modeler_name) in zip(modelers, modelers_names)
                for (solver, solver_name) in zip(solvers, solvers_names)
                    Test.@testset "$(solver_name) $(modeler_name)" verbose=VERBOSE showtiming=SHOWTIMING begin
                        sol = CommonSolve.solve(docp, init, modeler, solver; display=false)
                        Test.@test sol isa CTModels.Solution
                        Test.@test CTModels.successful(sol)
                        Test.@test CTModels.objective(sol) ≈ beam_data.obj atol=1e-2
                        Test.@test CTModels.iterations(sol) <= ipopt_options[:max_iter]
                        Test.@test CTModels.constraints_violation(sol) <= ipopt_options[:tol]
                    end 
                end
            end
        end=#

        # check_problem toplevel function
        Test.@testset "test_problem goddard" verbose=VERBOSE showtiming=SHOWTIMING begin
            test_problem(goddard(); display=true)
            test_problem(goddard2(); modeler=:exa, display=true) 
        end




    end

end