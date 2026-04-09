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

function build_exact_stagewise_xu(docp::CTDirect.DOCP{<:CTDirect.Gauss_Legendre_2_Stagewise})
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

function test_gauss_legendre_2_stagewise()
    prob = stagewise_scalar_problem()
    grid = [0.0, 0.2, 0.6, 1.0]

    @testset verbose = true showtiming = true ":gauss_legendre_2_stagewise :docp" begin
        docp = CTDirect.DOCP(prob.ocp, length(grid) - 1, 1, :gauss_legendre_2_stagewise, grid)
        disc = CTDirect.disc_model(docp)

        @test disc isa CTDirect.Gauss_Legendre_2_Stagewise
        @test docp.time.fixed_grid ≈ grid
        @test docp.time.control_steps == 1
        @test docp.dims.NLP_x == 1
        @test docp.dims.NLP_u == 1
        @test docp.dims.NLP_v == 0
        @test docp.dims.boundary_cons == 2
        @test disc.stage == 2
        @test disc._control_block == docp.dims.NLP_u
        @test disc._step_variables_block == 5
        @test disc._state_stage_eqs_block == 3
        @test disc._step_pathcons_block == 0
        @test docp.dim_NLP_variables == 16
        @test docp.dim_NLP_constraints == 11

        CTDirect.__variables_bounds!(docp)
        @test only(CTDirect.get_stagecontrol_at_time_step(docp.bounds.var_l, docp, 1, 1)) == 0.0
        @test only(CTDirect.get_stagecontrol_at_time_step(docp.bounds.var_u, docp, 1, 1)) == 2.0
        @test only(CTDirect.get_stagecontrol_at_time_step(docp.bounds.var_l, docp, 1, 2)) == 0.0
        @test only(CTDirect.get_stagecontrol_at_time_step(docp.bounds.var_u, docp, 1, 2)) == 2.0
    end

    @testset verbose = true showtiming = true ":gauss_legendre_2_stagewise :exact_feasible_trajectory" begin
        docp = CTDirect.DOCP(prob.ocp, length(grid) - 1, 1, :gauss_legendre_2_stagewise, grid)
        disc = CTDirect.disc_model(docp)
        xu = build_exact_stagewise_xu(docp)

        u11 = only(CTDirect.get_stagecontrol_at_time_step(xu, docp, 1, 1))
        u12 = only(CTDirect.get_stagecontrol_at_time_step(xu, docp, 1, 2))
        t22 = grid[2] + disc.butcher_c[2] * (grid[3] - grid[2])

        @test u11 != u12
        @test only(CTDirect.get_OCP_control_at_time_step(xu, docp, 1)) ≈ u11 atol = 1e-12
        @test only(CTDirect.get_OCP_state_at_time_step(xu, docp, 2)) ≈ grid[2]^2 atol = 1e-12
        @test only(CTDirect.get_stagevars_at_time_step(xu, docp, 2, 2)) ≈ 2 * t22 atol = 1e-12

        CTDirect.__constraints_bounds!(docp)
        c = zeros(docp.dim_NLP_constraints)
        CTDirect.__constraints!(c, xu, docp)

        @test c ≈ docp.bounds.con_l atol = 1e-12
        @test c ≈ docp.bounds.con_u atol = 1e-12
        @test CTDirect.__objective(xu, docp) ≈ (4 / 3) atol = 1e-12
    end

    @testset verbose = true showtiming = true ":gauss_legendre_2_stagewise :solve" begin
        sol = solve_problem(prob; scheme=:gauss_legendre_2_stagewise, grid_size=20, max_iter=200)
        T = time_grid(sol, :state)

        @test CTModels.successful(sol)
        @test objective(sol) ≈ prob.obj rtol = 1e-2
        @test first(T) ≈ 0.0 atol = 1e-12
        @test last(T) ≈ 1.0 atol = 1e-12
        @test only(state(sol)(first(T))) ≈ 0.0 atol = 1e-8
        @test only(state(sol)(last(T))) ≈ 1.0 atol = 1e-4
    end

    @testset verbose = true showtiming = true ":gauss_legendre_2_stagewise :compare_with_gauss_legendre_2" begin
        compare_grid = collect(LinRange(0.0, 1.0, 21))
        warmup_grid = collect(LinRange(0.0, 1.0, 5))
        docp_gl2 = CTDirect.DOCP(prob.ocp, length(compare_grid) - 1, 1, :gauss_legendre_2, compare_grid)
        docp_stagewise = CTDirect.DOCP(
            prob.ocp,
            length(compare_grid) - 1,
            1,
            :gauss_legendre_2_stagewise,
            compare_grid,
        )

        @test docp_gl2.time.fixed_grid ≈ docp_stagewise.time.fixed_grid
        @test docp_gl2.dim_NLP_constraints == docp_stagewise.dim_NLP_constraints
        @test docp_stagewise.dim_NLP_variables ==
            docp_gl2.dim_NLP_variables +
            docp_stagewise.time.steps * docp_stagewise.dims.NLP_u * docp_stagewise.time.control_steps

        solve_problem(prob; scheme=:gauss_legendre_2, time_grid=warmup_grid, max_iter=1)
        solve_problem(prob; scheme=:gauss_legendre_2_stagewise, time_grid=warmup_grid, max_iter=1)

        sol_gl2, elapsed_gl2 = solve_problem_timed(
            prob;
            scheme=:gauss_legendre_2,
            time_grid=compare_grid,
            max_iter=200,
        )
        sol_stagewise, elapsed_stagewise = solve_problem_timed(
            prob;
            scheme=:gauss_legendre_2_stagewise,
            time_grid=compare_grid,
            max_iter=200,
        )

        println(
            "gauss_legendre_2: vars=$(docp_gl2.dim_NLP_variables), cons=$(docp_gl2.dim_NLP_constraints), " *
            "iters=$(iterations(sol_gl2)), elapsed=$(round(elapsed_gl2; digits=3))s",
        )
        println(
            "gauss_legendre_2_stagewise: vars=$(docp_stagewise.dim_NLP_variables), cons=$(docp_stagewise.dim_NLP_constraints), " *
            "iters=$(iterations(sol_stagewise)), elapsed=$(round(elapsed_stagewise; digits=3))s",
        )

        @test CTModels.successful(sol_gl2)
        @test CTModels.successful(sol_stagewise)
        @test objective(sol_gl2) ≈ prob.obj rtol = 1e-2
        @test objective(sol_stagewise) ≈ prob.obj rtol = 1e-2
        @test objective(sol_stagewise) ≈ objective(sol_gl2) rtol = 1e-4 atol = 1e-6
        @test time_grid(sol_gl2, :state) ≈ time_grid(sol_stagewise, :state)
        @test only(state(sol_gl2)(1.0)) ≈ only(state(sol_stagewise)(1.0)) atol = 1e-4
    end
end
