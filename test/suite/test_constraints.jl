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
    ocp = goddard_a()
    sol = direct_solve(ocp.ocp, display=false, init=ocp.init)

    # plot state, control and costate
    psol = plot(sol)

    # check constraints and multipliers
    T = sol.times

    # POINT CONSTRAINTS
    # boundary constraints (NB. ORDERING IS NOT OBVIOUS)
    println("Boundary constraints: ", sol.infos[:boundary_constraints])
    println("multipliers: ", sol.infos[:mult_boundary_constraints])
    println("p(t0): ", sol.costate(T[1]), " p_m(tf) ", sol.costate(T[end])[3])
    # variable constraints / box
    println("\nVariable constraints: ", sol.infos[:variable_constraints])
    println("multipliers: ", sol.infos[:mult_variable_constraints])
    println("\nVariable box multipliers LB: ", sol.infos[:mult_variable_box_lower], " UB ", sol.infos[:mult_variable_box_upper])

    # PATH CONSTRAINTS
    # state box
    x = invert(sol.state.(T))
    x_box_lb = invert(sol.infos[:mult_state_box_lower].(T))
    x_box_ub = invert(sol.infos[:mult_state_box_upper].(T))

    p1_a = plot(T, x[1], label="state r")    
    p1_b = plot(T, [x_box_lb[1] x_box_ub[1]], label=["LB multiplier" "UB multiplier"])
    p1 = plot(p1_a, p1_b, layout=(2,1))

    p2_a = plot(T, x[2], label="state v")
    p2_b = plot(T, [x_box_lb[2] x_box_ub[2]], label=["LB multiplier" "UB multiplier"])
    p2 = plot(p2_a, p2_b, layout=(2,1))

    p3_a = plot(T, x[3], label="state m")
    p3_b = plot(T, [x_box_lb[3] x_box_ub[3]], label=["LB multiplier" "UB multiplier"])
    p3 = plot(p3_a, p3_b, layout=(2,1))

    # control box
    u = sol.control.(T)
    u_box_lb = flatten(sol.infos[:mult_control_box_lower].(T))
    u_box_ub = flatten(sol.infos[:mult_control_box_upper].(T))
    p4_a = plot(T, u, label="control u")    
    p4_b = plot(T, [u_box_lb u_box_ub], label=["LB multiplier" "UB multiplier"])
    p4 = plot(p4_a, p4_b, layout=(2,1))

    # display all graphs together
    p = plot(p1, p2, p3, p4, title=["r box" "" "v box" "" "m box" "" "u box" ""])
    display(p)

    # nonlinear path constraints
    # control constraints
    # state constraints
    # mixed constraints


end