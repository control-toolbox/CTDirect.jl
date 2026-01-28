println("testing: initial guess and continuation")

#################################################
# 1 Pass initial guess to all-in-one solve call

# use 0 iterations to check initial guess, >0 to check cv
maxiter = 0

# reference solution
if !isdefined(Main, :double_integrator_mintf)
    include("../problems/double_integrator.jl")
end
prob = double_integrator_mintf()
prob2 = double_integrator_minenergy()
sol0 = solve_problem(prob)

# constant initial guess
x_const = [0.5, 0.2]
u_const = 0.5
v_const = 0.15

# functional initial guess
x_func = t -> [t^2, sqrt(t)]
u_func = t -> (cos(10 * t) + 1) * 0.5

# interpolated initial guess
x_vec = [[0, 0], [1, 2], [5, -1]]
x_matrix = [0 0; 1 2; 5 -1]
u_vec = [0, 0.3, 0.1]

# 1.a default initial guess
@testset verbose = true showtiming = true ":default_init_no_arg" begin
    sol = solve_problem(prob; max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), (t -> [0.1, 0.1]).(T), rtol=1e-2)
    @test isapprox(control(sol).(T), (t -> 0.1).(T), rtol=1e-2)
    @test isapprox(variable(sol), 0.1, rtol=1e-2)
end
@testset verbose = true showtiming = true ":default_init_()" begin
    sol = solve_problem(prob; init=(), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), (t -> [0.1, 0.1]).(T), rtol=1e-2)
    @test isapprox(control(sol).(T), (t -> 0.1).(T), rtol=1e-2)
    @test isapprox(variable(sol), 0.1, rtol=1e-2)
end
# +++ NB this one will now use the initial guess provided in the problem itself...
# but this is the behaviour of local solve_problem, check vs OptimalControl solve !
@testset verbose = true showtiming = true ":default_init_nothing" begin
    sol = solve_problem(prob; init=nothing, max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), (t -> [0.1, 0.1]).(T), rtol=1e-2)
    @test isapprox(control(sol).(T), (t -> 0.1).(T), rtol=1e-2)
    @test variable(sol) == 0.1
end

# 1.b constant initial guess
@testset verbose = true showtiming = true ":constant_x" begin
    sol = solve_problem(prob; init=(state=x_const,), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), (t -> x_const).(T), rtol=1e-2)
end
@testset verbose = true showtiming = true ":constant_u" begin
    sol = solve_problem(prob; init=(control=u_const,), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(control(sol).(T), (t -> u_const).(T), rtol=1e-2)
end
@testset verbose = true showtiming = true ":constant_v" begin
    sol = solve_problem(prob; init=(variable=v_const,), max_iter=maxiter)
    @test variable(sol) == v_const
end
@testset verbose = true showtiming = true ":constant_xu" begin
    sol = solve_problem(prob; init=(state=x_const, control=u_const), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), (t -> x_const).(T), rtol=1e-2)
    @test isapprox(control(sol).(T), (t -> u_const).(T), rtol=1e-2)
end
@testset verbose = true showtiming = true ":constant_xu" begin
    sol = solve_problem(prob2; init=(state=x_const, control=u_const), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), (t -> x_const).(T), rtol=1e-2)
    @test isapprox(control(sol).(T), (t -> u_const).(T), rtol=1e-2)
end
@testset verbose = true showtiming = true ":constant_xv" begin
    sol = solve_problem(prob; init=(state=x_const, variable=v_const), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), (t -> x_const).(T), rtol=1e-2)
    @test variable(sol) == v_const
end
@testset verbose = true showtiming = true ":constant_uv" begin
    sol = solve_problem(prob; init=(control=u_const, variable=v_const), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(control(sol).(T), (t -> u_const).(T), rtol=1e-2)
    @test variable(sol) == v_const
end
@testset verbose = true showtiming = true ":constant_xuv" begin
    sol = solve_problem(prob;
        init=(state=x_const, control=u_const, variable=v_const),
        max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), (t -> x_const).(T), rtol=1e-2)
    @test isapprox(control(sol).(T), (t -> u_const).(T), rtol=1e-2)
    @test variable(sol) == v_const
end

# 1. functional initial guess
@testset verbose = true showtiming = true ":functional_x" begin
    sol = solve_problem(prob; init=(state=x_func,), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), x_func.(T), rtol=1e-2)
end
@testset verbose = true showtiming = true ":functional_u" begin
    sol = solve_problem(prob; init=(control=u_func,), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(control(sol).(T), u_func.(T), rtol=1e-2)
end
@testset verbose = true showtiming = true ":functional_xu" begin
    sol = solve_problem(prob; init=(state=x_func, control=u_func), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), x_func.(T), rtol=1e-2)
    @test isapprox(control(sol).(T), u_func.(T), rtol=1e-2)
end
@testset verbose = true showtiming = true ":functional_xu" begin
    sol = solve_problem(prob2; init=(state=x_func, control=u_func), max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), x_func.(T), rtol=1e-2)
    @test isapprox(control(sol).(T), u_func.(T), rtol=1e-2)
end

# 1.d interpolated initial guess
t_vec = [0, 0.1, v_const]
@testset verbose = true showtiming = true ":vector_txu :constant_v" begin
    init = (state=(t_vec, x_vec), control=(t_vec, u_vec), variable=v_const)
    sol = solve_problem(prob; init=init, max_iter=maxiter)
    @test isapprox(state(sol).(t_vec), x_vec, rtol=1e-2)
    @test isapprox(control(sol).(t_vec), u_vec, rtol=1e-2)
    @test variable(sol) == v_const
end
t_matrix = [0 0.1 v_const]
@testset verbose = true showtiming = true ":matrix_t :vector_xu :constant_v" begin
    init = (state=(t_matrix, x_vec), control=(t_matrix, u_vec), variable=v_const)
    sol = solve_problem(prob; init=init, max_iter=maxiter)
    @test isapprox(state(sol).(flatten(t_matrix)), x_vec, rtol=1e-2)
    @test isapprox(control(sol).(flatten(t_matrix)), u_vec, rtol=1e-2)
    @test variable(sol) == v_const
end
@testset verbose = true showtiming = true ":matrix_x :vector_tu :constant_v" begin
    init = (state=(t_vec, x_matrix), control=(t_vec, u_vec), variable=v_const)
    sol = solve_problem(prob; init=init, max_iter=maxiter)
    @test isapprox(stack(state(sol).(t_matrix), dims=1), x_matrix, rtol=1e-2)
    @test isapprox(control(sol).(t_vec), u_vec, rtol=1e-2)
    @test variable(sol) == v_const
end

# 1.e mixed initial guess
@testset verbose = true showtiming = true ":vector_tx :functional_u :constant_v" begin
    init = (state=(t_vec, x_vec), control=u_func, variable=v_const)
    sol = solve_problem(prob; init=init, max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(t_vec), x_vec, rtol=1e-2)
    @test isapprox(control(sol).(T), u_func.(T), rtol=1e-2)
    @test variable(sol) == v_const
end

# 1.f warm start
@testset verbose = true showtiming = true ":warm_start" begin
    sol = solve_problem(prob; init=sol0, max_iter=maxiter)
    T = time_grid(sol)
    @test isapprox(state(sol).(T), state(sol0).(T), rtol=1e-2)
    @test isapprox(control(sol).(T), control(sol0).(T), rtol=1e-2)
    @test variable(sol) == variable(sol0)
end


##########################
# 2. discrete continuation

test1 = true
test2 = true
test3 = true
draw_plot = false

# 2.1 double integrator 
if test1
    if !isdefined(Main, :double_integrator_minenergy)
        include("../problems/double_integrator.jl")
    end
    @testset verbose = true showtiming = true ":continuation :double_integrator" begin
        init = ()
        obj_list = []
        for T in 1:5
            prob = double_integrator_minenergy(T)
            sol = solve_problem(prob; init=init, grid=100)
            init = sol
            push!(obj_list, objective(sol))
        end
        @test obj_list ≈ [12, 1.5, 0.44, 0.19, 0.096] rtol = 1e-2
    end
end

# 2.2 parametric
if test2
    if !isdefined(Main, :parametric)
        include("../problems/parametric.jl")
    end
    @testset verbose = true showtiming = true ":continuation :parametric_ocp" begin
        init = ()
        obj_list = []
        for ρ in [0.1, 5, 10, 30, 100]
            prob = parametric(ρ)
            sol = solve_problem(prob; init=init)
            init = sol
            push!(obj_list, objective(sol))
        end
        @test obj_list ≈ [-0.034, -1.67, -6.2, -35, -148] rtol = 1e-2
    end
end

# 2.3 goddard
if test3
    if !isdefined(Main, :goddard)
        include("../problems/goddard.jl")
    end
    sol0 = solve_problem(goddard())
    @testset verbose = true showtiming = true ":continuation :goddard" begin
        sol = sol0
        Tmax_list = []
        obj_list = []
        for Tmax in 3.5:-0.5:1
            sol = solve_problem(goddard(Tmax=Tmax); init=sol)
            push!(Tmax_list, Tmax)
            push!(obj_list, objective(sol))
        end
        @test obj_list ≈ [1.0125, 1.0124, 1.0120, 1.0112, 1.0092, 1.0036] rtol = 1e-2

        if draw_plot
            using Plots
            # plot obj(vmax)
            pobj = plot(
                Tmax_list,
                obj_list,
                label="r(tf)",
                xlabel="Maximal thrust (Tmax)",
                ylabel="Maximal altitude r(tf)",
                seriestype=:scatter,
            )
            # plot multiple solutions
            plot(sol0)
            p = plot!(sol)
            display(plot(pobj, p, layout=2, reuse=false, size=(1000, 500)))
        end
    end
end
