# Direct methods

"""
$(TYPEDSIGNATURES)

Default discretization method: `midpoint`.
"""
__disc_method() = :midpoint

"""
$(TYPEDSIGNATURES)

Default grid size: `250`.
"""
__grid_size() = 250

"""
$(TYPEDSIGNATURES)

Default initial guess: `nothing` (corresponds to 0.1 for all variables)
"""
__ocp_init() = nothing

"""
$(TYPEDSIGNATURES)

Default display toggle: `true`
"""
__display() = true

"""
$(TYPEDSIGNATURES)

Default (non uniform) time grid: `nothing`
"""
__time_grid() = nothing

"""
$(TYPEDSIGNATURES)

Default backend for ADNLPModels: `:optimized`
"""
__adnlp_backend() = :optimized

"""
$(TYPEDSIGNATURES)

Default backend for ExaModels: `nothing`
"""
__exa_backend() = nothing

"""
$(TYPEDSIGNATURES)

Reformulate Lagrange cost as Mayer cost: false
"""
__lagrange_to_mayer() = false

"""
$(TYPEDSIGNATURES)

Default tolerance: `1e-8`
"""
__tolerance() = 1e-8

"""
$(TYPEDSIGNATURES)

Default maximum of iterations: `1000`
"""
__max_iterations() = 1000


