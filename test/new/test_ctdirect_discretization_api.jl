# Unit tests for the discretization API (discretize with custom and default discretizers).
struct DummyOCPDiscretize <: CTModels.AbstractOptimalControlProblem end

struct DummyDiscretizer <: CTDirect.AbstractOptimalControlDiscretizer
    calls::Base.RefValue{Int}
    tag::Symbol
end

function (d::DummyDiscretizer)(ocp::CTModels.AbstractOptimalControlProblem)
    d.calls[] += 1
    return (ocp, d.tag)
end

function test_ctdirect_discretization_api()

    # ========================================================================
    # discretize(ocp, discretizer)
    # ========================================================================

    Test.@testset "discretize(ocp, discretizer)" verbose=VERBOSE showtiming=SHOWTIMING begin
        ocp = DummyOCPDiscretize()
        calls = Ref(0)
        discretizer = DummyDiscretizer(calls, :dummy)

        result = CTDirect.discretize(ocp, discretizer)

        Test.@test result == (ocp, :dummy)
        Test.@test calls[] == 1
    end

    # ========================================================================
    # discretize(ocp; discretizer=__discretizer())
    # ========================================================================

    Test.@testset "default discretizer" verbose=VERBOSE showtiming=SHOWTIMING begin
        ocp = DummyOCPDiscretize()

        docp = CTDirect.discretize(ocp)

        # The default discretizer should produce a DiscretizedOptimalControlProblem
        Test.@test docp isa CTModels.DiscretizedOptimalControlProblem
        Test.@test CTModels.ocp_model(docp) === ocp

        # And the low-level __discretizer() helper should return a Collocation
        disc = CTDirect.__discretizer()
        Test.@test disc isa CTDirect.AbstractOptimalControlDiscretizer
        Test.@test disc isa CTDirect.Collocation
    end
end