println("testing: path and boundary duals by label")

# labelled path constraints (one and two rows), labelled boundary constraints, a control box
function labelled_constraints()
    @def ocp begin
        t ∈ [0, 1], time
        x ∈ R², state
        u ∈ R, control
        -2 ≤ u(t) ≤ 2, (u_box)
        x(0) == [-1, 0], (x_init)
        x₁(1) + x₂(1)^2 == 0, (x_final)
        x₂(t) + u(t) ≤ 1.2, (c_one)
        [x₂(t) - u(t), x₁(t)^2] ≤ [2, 4], (c_two)
        ẋ(t) == [x₂(t), u(t)]
        ∫(u(t)^2) → min
    end

    return ((ocp=ocp, obj=nothing, name="labelled_constraints", init=()))
end

const PATH_LABELS = (:c_one, :c_two)
const BOUNDARY_LABELS = (:x_init, :x_final)

@testset verbose = true showtiming = true ":duals_by_label :adnlp" begin
    prob = labelled_constraints()
    ocp = prob.ocp
    sol = @test_nowarn solve_problem(prob; grid_size=50)
    @test CTModels.successful(sol)

    # the active constraint has a nonzero dual: the check below compares real values
    @test maximum(abs, CTModels.Solutions.dual(sol, ocp, :c_one).(time_grid(sol, :path))) > 1e-3

    pc = CTModels.path_constraints_dual(sol)
    for label in PATH_LABELS
        @test CTModels.Solutions.has_dual(sol, ocp, label)
        _, rows = CTModels.Models.constraint_rows(ocp, label)
        d = CTModels.Solutions.dual(sol, ocp, label)
        for t in time_grid(sol, :path)
            @test d(t) ≈ (length(rows) == 1 ? pc(t)[only(rows)] : pc(t)[rows])
        end
    end

    bc = CTModels.boundary_constraints_dual(sol)
    for label in BOUNDARY_LABELS
        @test CTModels.Solutions.has_dual(sol, ocp, label)
        _, rows = CTModels.Models.constraint_rows(ocp, label)
        expected = length(rows) == 1 ? bc[only(rows)] : bc[rows]
        @test CTModels.Solutions.dual(sol, ocp, label) ≈ expected
        @test all(isfinite, expected)
    end
end

@testset verbose = true showtiming = true ":duals_by_label :exa" begin
    prob = labelled_constraints()
    ocp = prob.ocp
    sol = @test_nowarn solve_problem(prob; modeler=:exa, grid_size=50)
    @test CTModels.successful(sol)

    # not returned by the exa getter yet: unknown, never 0
    for label in (PATH_LABELS..., BOUNDARY_LABELS...)
        @test !CTModels.Solutions.has_dual(sol, ocp, label)
        @test_throws Exceptions.PreconditionError CTModels.Solutions.dual(sol, ocp, label)
    end
    err = try
        CTModels.Solutions.dual(sol, ocp, :c_one)
    catch e
        e
    end
    @test occursin("CTParser.jl/issues/414", sprint(showerror, err))

    # the box duals are still given
    @test CTModels.Solutions.has_dual(sol, ocp, :u_box)

    # plots: the path constraints are drawn without their duals
    @test !isnothing(CTModels.PlotCase.build_figure(sol, :path, :dual))
end
