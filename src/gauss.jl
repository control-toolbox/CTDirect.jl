#= Functions for implicit gauss-legendre 2 discretization scheme
Internal layout for NLP variables: 
[X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the conventions U_i constant per step, K_i := [K_i1 K_i2]

NB. try both control_stage and control_step versions !
Q. if taking controls at stages, how to define the control at 'step' (for path constraints or tf). Unless we evaluate path constraints at stages, but then we have to determine the state at stages (a possibility is to take the state argument as for the stage dynamics)
=#

# Later adjust this one for generic IRK scheme ! (use stage value)
# define several tags for different methods
# use intermediate abstract type ImplicitRungeKuttaTag

struct GaussLegendre2Tag <: DiscretizationTag
    stage::Int
    additional_controls::Int
    butcher_a
    butcher_b
    butcher_c
    GaussLegendre2Tag() = new(
        2,
        0,
        [0.25 (0.25-sqrt(3) / 6); (0.25+sqrt(3) / 6) 0.25],
        [0.5, 0.5],
        [(0.5 - sqrt(3) / 6), (0.5 + sqrt(3) / 6)],
    )
end

"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
"""
function get_variables_at_time_step(xu, docp, i, tag::GuassLegendre2Tag)

    # block: [X_i U_i1 U_i2 K_i1 K_i2]
    nx = docp.dim_NLP_x
    n = docp.dim_OCP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx * 3 + m) * i

    # retrieve scalar/vector OCP state (w/o lagrange state) 
    if n == 1
        xi = xu[offset + 1]
    else
        xi = xu[(offset + 1):(offset + n)]
    end
    if docp.has_lagrange
        xli = xu[offset + nx]
    else
        xli = nothing # dummy. use xu type ?
    end

    # retrieve scalar/vector controls
    if i < N
        offset_u = offset
    else
        offset_u = (nx * 3 + m) * (i - 1)
    end
    if m == 1
        ui = xu[offset_u + nx + 1]
    else
        ui = xu[(offset_u + nx + 1):(offset_u + nx + m)]
    end

    # retrieve vector stage variable (except at final time)
    if i < N
        ki = (
            xu[(offset + nx + m + 1):(offset + nx + m + nx)],
            xu[(offset + nx * 2 + m + 1):(offset + nx * 2 + m + nx)],
        )
    else
        ki = nothing
    end

    return xi, ui, xli, ki
end

# internal NLP version for solution parsing
# could be fused with one above if 
# - using extended dynamics that include lagrange cost
# - scalar case is handled at OCP level
function get_NLP_variables_at_time_step(xu, docp, i, tag::GaussLegendre2Tag)
    ++ + nx = docp.dim_NLP_x
    m = docp.dim_NLP_u
    N = docp.dim_NLP_steps
    offset = (nx * 2 + m) * i

    # state
    xi = xu[(offset + 1):(offset + nx)]
    # control
    if i < N
        offset_u = offset
    else
        offset_u = (nx * 2 + m) * (i - 1)
    end
    ui = xu[(offset_u + nx + 1):(offset_u + nx + m)]
    # stage
    if i < N
        ki = xu[(offset + nx + m + 1):(offset + nx + m + nx)]
    else
        ki = nothing
    end

    return xi, ui, ki
end
