println("Test: constraint types")

if !isdefined(Main, :goddard)
    include("../problems/goddard.jl")
end

# constraints / multipliers check
check_constraint_mult = true 


#=
# box constraints
@testset verbose = true showtiming = true ":goddard :box_constraints" begin
    ocp = goddard()
    sol = direct_solve(ocp.ocp, display=false)
    @test sol.objective ≈ ocp.obj rtol=1e-2
end

# functional constraints
@testset verbose = true showtiming = true ":goddard :functional_constraints" begin
    ocp = goddard(functional_constraints=true)
    sol = direct_solve(ocp.ocp, display=false, init=ocp.init)
    @test sol.objective ≈ ocp.obj rtol=1e-2
end

# all constraints
@testset verbose = true showtiming = true ":goddard :all_constraints" begin
    ocp = goddard_all()
    sol = direct_solve(ocp.ocp, display=false, init=ocp.init)
    @test sol.objective ≈ ocp.obj rtol=1e-2
end
=#


# NB. sign mismatch between p(tf) and multiplier for final constraints 
# Note: if possible fix ordering in OCP: initial THEN final constraints !

if check_constraint_mult
    using SplitApplyCombine
    ocp = goddard_all()
    sol = direct_solve(ocp.ocp, display=false, init=ocp.init)

    # plot state, control and costate
    psol = plot(sol)

    # check constraints and multipliers
    T = sol.times
    
    # +++ parse and plot constraints

    # boundary constraints (NB. ORDERING IS NOT OBVIOUS)
    println("Boundary constraints: ", sol.infos[:boundary_constraints])
    println("multipliers: ", sol.infos[:mult_boundary_constraints])
    println("p(t0): ", sol.costate(T[1]), " p_m(tf) ", sol.costate(T[end])[3])
    # variable constraints / box
    println("\nVariable constraints: ", sol.infos[:variable_constraints])
    println("multipliers: ", sol.infos[:mult_variable_constraints])
    println("\nVariable box multipliers LB: ", sol.infos[:mult_variable_box_lower], " UB ", sol.infos[:mult_variable_box_upper])

    # box constraints
    # state box
    # NB SCALE IS BAD plot variable then multipliers below (split 2x1 graph)
    x = invert(sol.state.(T))
    x_box_lb = invert(sol.infos[:mult_state_box_lower].(T))
    x_box_ub = invert(sol.infos[:mult_state_box_upper].(T))
    p1 = plot(T, x[1], label="state r")    
    p1 = plot!(T, x_box_lb[1], label="LB multiplier")
    p1 = plot!(T, x_box_ub[1], label="UB multiplier")
    p2 = plot(T, x[2], label="state v")
    p2 = plot!(T, x_box_lb[2], label="LB multiplier")
    p2 = plot!(T, x_box_ub[2], label="UB multiplier")
    p3 = plot(T, x[3], label="state m")    
    p3 = plot!(T, x_box_lb[3], label="LB multiplier")
    p3 = plot!(T, x_box_ub[3], label="UB multiplier")
    display(p3)
    # control box
    u = sol.control.(T)
    u_box_lb = flatten(sol.infos[:mult_control_box_lower].(T))
    u_box_ub = flatten(sol.infos[:mult_control_box_upper].(T))
    p4 = plot(T, u, label="control u")    
    p4 = plot!(T, u_box_lb, label="LB multiplier")
    p4 = plot!(T, u_box_ub, label="UB multiplier")
    display(p4)

    # overall use 2x2 plot for x,u boxes


    # nonlinear path constraints


end