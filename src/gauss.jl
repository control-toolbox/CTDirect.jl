#= Functions for implicit gauss-legendre 2 discretization scheme
Internal layout for NLP variables: 
[X_0,U_0,K_0, X_1,U_1,K_1 .., X_N-1,U_N-1,K_N-1, XN, V]
with the conventions U_i constant per step, K_i := [K_i1 K_i2]

NB. try both control_stage and control_step versions !
Q. if taking controls at stages, how to define the control at 'step' (for path constraints or tf). Unless we evaluate path constraints at stages, but then we have to determine the state at stages (a possibility is to take the state argument as for the stage dynamics)
=#

# Later adjust this one for generic IRK scheme ! (use stage value)
# define several tags for different methods
# use intermediate abstract type ImplicitRungeKutta ?
#+++ add Midpoint_imp

struct GaussLegendre2 <: Discretization
    stage::Int
    additional_controls::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    GaussLegendre2() = new(
        2,
        0,
        [0.25 (0.25-sqrt(3) / 6); (0.25+sqrt(3) / 6) 0.25],
        [0.5, 0.5],
        [(0.5 - sqrt(3) / 6), (0.5 + sqrt(3) / 6)],
    )
end
