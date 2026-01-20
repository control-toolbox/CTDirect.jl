# Unit tests for Collocation discretizer wiring from OCP to discretized OCP and builders.
struct DummyOCPCollocation <: CTModels.AbstractOptimalControlProblem end

struct DummyOCPExaRouting <: CTModels.AbstractOptimalControlProblem end

struct DummyDOCPCollocationRouting end

const CM_ExaRecordedCollocation = Ref{Any}(nothing)

#=
function CTDirect.direct_transcription(
    ocp::DummyOCPExaRouting,
    modeler::Symbol;
    grid_size,
    disc_method,
    init,
    lagrange_to_mayer,
    kwargs...,
)
    CM_ExaRecordedCollocation[] = (
        ocp=ocp,
        modeler=modeler,
        grid_size=grid_size,
        disc_method=disc_method,
        init=init,
        lagrange_to_mayer=lagrange_to_mayer,
        kwargs=NamedTuple(kwargs),
    )
    return DummyDOCPCollocationRouting()
end

function CTDirect.nlp_model(::DummyDOCPCollocationRouting)
    # Build a minimal but well-formed ExaModels.ExaModel: one variable and a
    # trivial objective, no constraints. This exercises the collocation path
    # end-to-end without relying on a specific test problem.
    BaseType = Float64
    core = ExaModels.ExaCore(BaseType)
    x = ExaModels.variable(core, 1; start=BaseType[0])
    ExaModels.objective(core, x[1])
    return ExaModels.ExaModel(core)
end
=#

function test_ctdirect_collocation()
    Test.@testset "Collocation as discretizer" verbose=VERBOSE showtiming=SHOWTIMING begin
        ocp = DummyOCPCollocation()

        # Use the default Collocation discretizer to avoid relying on CTDirect
        discretizer = CTDirect.__discretizer()
        Test.@test discretizer isa CTDirect.Collocation

        docp = discretizer(ocp)

        # The call operator on Collocation should return a DiscretizedOptimalControlProblem
        Test.@test docp isa CTModels.DiscretizedOptimalControlProblem
        Test.@test CTModels.ocp_model(docp) === ocp

        # The model and solution builders should be correctly wired with both
        # ADNLP and Exa backends present.
        adnlp_builder = CTModels.get_adnlp_model_builder(docp)
        exa_builder = CTModels.get_exa_model_builder(docp)
        adnlp_sol = CTModels.get_adnlp_solution_builder(docp)
        exa_sol = CTModels.get_exa_solution_builder(docp)

        Test.@test adnlp_builder isa CTModels.ADNLPModelBuilder
        Test.@test exa_builder isa CTModels.ExaModelBuilder
        Test.@test adnlp_sol isa CTModels.ADNLPSolutionBuilder
        Test.@test exa_sol isa CTModels.ExaSolutionBuilder
    end

    Test.@testset "Exa backend routing" verbose=VERBOSE showtiming=SHOWTIMING begin
        ocp = DummyOCPExaRouting()

        # Stub CTDirect.direct_transcription for DummyOCPExaRouting to record kwargs
        CM_ExaRecordedCollocation[] = nothing

        # Case 1: default grid (Int) and default lagrange_to_mayer=false
        discretizer = CTDirect.Collocation()
        docp = discretizer(ocp)

        exa_builder = CTModels.get_exa_model_builder(docp)

        # Minimal initial guess: functions for state/control and empty variable
        init_guess = CTModels.OptimalControlInitialGuess(t -> 0.0, t -> 0.0, Float64[])

        BaseType = Float32
        exa_nlp = exa_builder(BaseType, init_guess; backend=:gpu, foo=1)
        Test.@test exa_nlp isa ExaModels.ExaModel

        # The direct_transcription stub must have recorded the call.
        rec = CM_ExaRecordedCollocation[]
        Test.@test rec !== nothing
        Test.@test rec[:modeler] === :exa

        grid_default = CTModels.get_option_value(discretizer, :grid)
        Test.@test rec[:grid_size] == grid_default
        Test.@test rec[:lagrange_to_mayer] === false

        kw = rec[:kwargs]
        # time_grid should be absent or nothing for Int grid
        if haskey(kw, :time_grid)
            Test.@test kw[:time_grid] === nothing
        end

        # backend should have been rerouted to exa_backend
        Test.@test haskey(kw, :exa_backend)
        Test.@test kw[:exa_backend] === :gpu
        # original backend key should not be forwarded
        Test.@test !haskey(kw, :backend)
        # other kwargs are preserved
        Test.@test haskey(kw, :foo)
        Test.@test kw[:foo] == 1

        # Case 2: explicit time grid (Vector) and lagrange_to_mayer=true
        CM_ExaRecordedCollocation[] = nothing
        grid_vec = collect(range(0.0, 1.0; length=5))
        discretizer2 = CTDirect.Collocation(; grid=grid_vec, lagrange_to_mayer=true)
        docp2 = discretizer2(ocp)
        exa_builder2 = CTModels.get_exa_model_builder(docp2)
        exa_nlp2 = exa_builder2(BaseType, init_guess; backend=:gpu)
        Test.@test exa_nlp2 isa ExaModels.ExaModel

        rec2 = CM_ExaRecordedCollocation[]
        Test.@test rec2 !== nothing
        Test.@test rec2[:modeler] === :exa
        Test.@test rec2[:grid_size] == length(grid_vec)
        Test.@test rec2[:lagrange_to_mayer] === true

        kw2 = rec2[:kwargs]
        Test.@test haskey(kw2, :time_grid)
        Test.@test kw2[:time_grid] === grid_vec
    end
end