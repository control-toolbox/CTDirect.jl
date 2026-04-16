function stagewise_scalar_problem()
    @def ocp begin
        t ∈ [0, 1], time
        x ∈ R, state
        u ∈ R, control
        0 ≤ u(t) ≤ 2
        x(0) == 0
        x(1) == 1
        ẋ(t) == u(t)
        ∫(u(t)^2) → min
    end

    return (ocp=ocp, obj=1.0, name="stagewise_scalar", init=())
end

stagewise_feasible_state(t) = [t^2]
stagewise_feasible_control(t) = [2 * t]

function build_exact_stagewise_xu(docp::CTDirect.DOCP{<:CTDirect.GenericIRKStagewise})
    xu = zeros(docp.dim_NLP_variables)
    disc = CTDirect.disc_model(docp)
    T = docp.time.fixed_grid

    for i in 1:(docp.time.steps + 1)
        CTDirect.get_OCP_state_at_time_step(xu, docp, i) .= stagewise_feasible_state(T[i])
    end

    for i in 1:docp.time.steps
        ti = T[i]
        hi = T[i + 1] - ti

        for j in 1:disc.stage
            tij = ti + disc.butcher_c[j] * hi
            CTDirect.get_stagecontrol_at_time_step(xu, docp, i, j) .=
                stagewise_feasible_control(tij)
            CTDirect.get_stagevars_at_time_step(xu, docp, i, j) .=
                stagewise_feasible_control(tij)
        end
    end

    return xu
end

function solve_problem_timed(prob; kwargs...)
    elapsed = @elapsed sol = solve_problem(prob; kwargs...)
    return sol, elapsed
end

function test_gauss_legendre_stagewise_scheme(spec)
    prob = stagewise_scalar_problem()
    grid = [0.0, 0.2, 0.6, 1.0]

    @testset verbose = true showtiming = true "$(spec.scheme) :docp" begin
        docp = CTDirect.DOCP(prob.ocp, length(grid) - 1, 1, spec.scheme, grid)
        disc = CTDirect.disc_model(docp)

        @test disc isa spec.disc_type
        @test docp.time.fixed_grid ≈ grid
        @test docp.time.control_steps == 1
        @test docp.dims.NLP_x == 1
        @test docp.dims.NLP_u == 1
        @test docp.dims.NLP_v == 0
        @test docp.dims.boundary_cons == 2
        @test disc.stage == spec.stage
        @test disc._control_block == docp.dims.NLP_u * spec.stage
        @test disc._step_variables_block == spec.expected_step_variables_block
        @test disc._state_stage_eqs_block == spec.expected_state_stage_eqs_block
        @test disc._step_pathcons_block == 0
        @test docp.dim_NLP_variables == spec.expected_dim_NLP_variables
        @test docp.dim_NLP_constraints == spec.expected_dim_NLP_constraints

        CTDirect.__variables_bounds!(docp)
        for j in 1:disc.stage
            @test only(CTDirect.get_stagecontrol_at_time_step(docp.bounds.var_l, docp, 1, j)) == 0.0
            @test only(CTDirect.get_stagecontrol_at_time_step(docp.bounds.var_u, docp, 1, j)) == 2.0
        end
    end

    @testset verbose = true showtiming = true "$(spec.scheme) :exact_feasible_trajectory" begin
        docp = CTDirect.DOCP(prob.ocp, length(grid) - 1, 1, spec.scheme, grid)
        disc = CTDirect.disc_model(docp)
        xu = build_exact_stagewise_xu(docp)

        u11 = only(CTDirect.get_stagecontrol_at_time_step(xu, docp, 1, 1))
        u1s = only(CTDirect.get_stagecontrol_at_time_step(xu, docp, 1, disc.stage))
        t2s = grid[2] + disc.butcher_c[disc.stage] * (grid[3] - grid[2])

        @test u11 != u1s
        @test only(CTDirect.get_OCP_control_at_time_step(xu, docp, 1)) ≈ u11 atol = 1e-12
        @test only(CTDirect.get_OCP_state_at_time_step(xu, docp, 2)) ≈ grid[2]^2 atol = 1e-12
        @test only(CTDirect.get_stagevars_at_time_step(xu, docp, 2, disc.stage)) ≈ 2 * t2s atol = 1e-12

        CTDirect.__constraints_bounds!(docp)
        c = zeros(docp.dim_NLP_constraints)
        CTDirect.__constraints!(c, xu, docp)

        @test c ≈ docp.bounds.con_l atol = 1e-12
        @test c ≈ docp.bounds.con_u atol = 1e-12
        @test CTDirect.__objective(xu, docp) ≈ (4 / 3) atol = 1e-12
    end

    @testset verbose = true showtiming = true "$(spec.scheme) :solve" begin
        sol = solve_problem(prob; scheme=spec.scheme, grid_size=20, max_iter=200)
        T = time_grid(sol, :state)

        @test CTModels.successful(sol)
        @test objective(sol) ≈ prob.obj rtol = 1e-2
        @test first(T) ≈ 0.0 atol = 1e-12
        @test last(T) ≈ 1.0 atol = 1e-12
        @test only(state(sol)(first(T))) ≈ 0.0 atol = 1e-8
        @test only(state(sol)(last(T))) ≈ 1.0 atol = 1e-4
    end

    @testset verbose = true showtiming = true "$(spec.scheme) :compare_with_$(spec.reference_scheme)" begin
        compare_grid = collect(LinRange(0.0, 1.0, 21))
        warmup_grid = collect(LinRange(0.0, 1.0, 5))
        docp_reference = CTDirect.DOCP(
            prob.ocp,
            length(compare_grid) - 1,
            1,
            spec.reference_scheme,
            compare_grid,
        )
        docp_stagewise = CTDirect.DOCP(
            prob.ocp,
            length(compare_grid) - 1,
            1,
            spec.scheme,
            compare_grid,
        )

        @test docp_reference.time.fixed_grid ≈ docp_stagewise.time.fixed_grid
        @test docp_reference.dim_NLP_constraints == docp_stagewise.dim_NLP_constraints
        @test docp_stagewise.dim_NLP_variables ==
            docp_reference.dim_NLP_variables +
            docp_stagewise.time.steps *
            (CTDirect.disc_model(docp_stagewise).stage - 1) *
            docp_stagewise.dims.NLP_u *
            docp_stagewise.time.control_steps

        solve_problem(prob; scheme=spec.reference_scheme, time_grid=warmup_grid, max_iter=1)
        solve_problem(prob; scheme=spec.scheme, time_grid=warmup_grid, max_iter=1)

        sol_reference, elapsed_reference = solve_problem_timed(
            prob;
            scheme=spec.reference_scheme,
            time_grid=compare_grid,
            max_iter=200,
        )
        sol_stagewise, elapsed_stagewise = solve_problem_timed(
            prob;
            scheme=spec.scheme,
            time_grid=compare_grid,
            max_iter=200,
        )

        println(
            "$(spec.reference_scheme): vars=$(docp_reference.dim_NLP_variables), cons=$(docp_reference.dim_NLP_constraints), " *
            "iters=$(iterations(sol_reference)), elapsed=$(round(elapsed_reference; digits=3))s",
        )
        println(
            "$(spec.scheme): vars=$(docp_stagewise.dim_NLP_variables), cons=$(docp_stagewise.dim_NLP_constraints), " *
            "iters=$(iterations(sol_stagewise)), elapsed=$(round(elapsed_stagewise; digits=3))s",
        )

        @test CTModels.successful(sol_reference)
        @test CTModels.successful(sol_stagewise)
        @test objective(sol_reference) ≈ prob.obj rtol = 1e-2
        @test objective(sol_stagewise) ≈ prob.obj rtol = 1e-2
        @test objective(sol_stagewise) ≈ objective(sol_reference) rtol = 1e-4 atol = 1e-6
        @test time_grid(sol_reference, :state) ≈ time_grid(sol_stagewise, :state)
        @test only(state(sol_reference)(1.0)) ≈ only(state(sol_stagewise)(1.0)) atol = 1e-4
    end
end

function test_stagewise()
    specs = (
        (
            scheme=:gauss_legendre_2,
            reference_scheme=:gauss_legendre_2_constant_control,
            disc_type=CTDirect.Gauss_Legendre_2_Stagewise,
            stage=2,
            expected_step_variables_block=5,
            expected_state_stage_eqs_block=3,
            expected_dim_NLP_variables=16,
            expected_dim_NLP_constraints=11,
        ),
        (
            scheme=:gauss_legendre_3,
            reference_scheme=:gauss_legendre_3_constant_control,
            disc_type=CTDirect.Gauss_Legendre_3_Stagewise,
            stage=3,
            expected_step_variables_block=7,
            expected_state_stage_eqs_block=4,
            expected_dim_NLP_variables=22,
            expected_dim_NLP_constraints=14,
        ),
    )

    for spec in specs
        test_gauss_legendre_stagewise_scheme(spec)
    end
end

test_stagewise()

