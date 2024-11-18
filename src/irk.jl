#= Functions for generic implicit Runge Kutta discretization
Internal layout for NLP variables: 
[X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the conventions U_i constant per step, K_i := [K_i1 K_i2]
=#

abstract type GenericIRK <: Discretization end

struct Midpoint_IRK <: GenericIRK
    stage::Int
    additional_controls::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    function Midpoint_IRK() = return new(
        1, 
        0, 
        hcat(0.5), 
        [1], 
        [0.5], 
        "Implicit Midpoint aka Gauss-Legendre collocation for s=1, 2nd order, symplectic")
end


struct GaussLegendre2 <: GenericIRK
    stage::Int
    additional_controls::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    function GaussLegendre2() = return new(
        2,
        0,
        [0.25 (0.25-sqrt(3) / 6); (0.25+sqrt(3) / 6) 0.25],
        [0.5, 0.5],
        [(0.5 - sqrt(3) / 6), (0.5 + sqrt(3) / 6)],
        "Implicit Gauss-Legendre collocation for s=2, 4th order, symplectic")
    )
end


"""
$(TYPEDSIGNATURES)

Retrieve state and control variables at given time step from the NLP variables.
Convention: 1 <= i <= dim_NLP_steps+1
"""
function get_OCP_state_at_time_step(xu, docp::DOCP{GenericIRK, ScalVariable, <: ScalVect, <: ScalVect}, i)
end
function get_OCP_state_at_time_step(xu, docp::DOCP{GenericIRK, VectVariable, <: ScalVect, <: ScalVect}, i)
end
function get_lagrange_state_at_time_step(xu, docp::DOCP{GenericIRK}, i)
end

function get_OCP_control_at_time_step(xu, docp::DOCP{GenericIRK, <: ScalVect, ScalVariable, <: ScalVect}, i)
end
function get_OCP_control_at_time_step(xu, docp::DOCP{GenericIRK, <: ScalVect, VectVariable, <: ScalVect}, i)
end

function get_stagevars_at_time_step(xu, docp::DOCP{GenericIRK}, i)
end


function set_state_at_time_step!(xu, x_init, docp::DOCP{GenericIRK}, i)
end
function set_control_at_time_step!(xu, u_init, docp::DOCP{GenericIRK}, i)
end


function setWorkArray(docp::DOCP{Midpoint}, xu, time_grid, v)
end

"""
$(TYPEDSIGNATURES)

Set the constraints corresponding to the state equation
Convention: 1 <= i <= dim_NLP_steps (+1)
"""
function setConstraintBlock!(docp::DOCP{Midpoint}, c, xu, v, time_grid, i, work)
end