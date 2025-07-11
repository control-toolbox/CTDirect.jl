# Direct methods

"""
$(TYPEDSIGNATURES)

Default discretization method: `trapeze`.
"""
__disc_method() = :trapeze

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

Default tolerance: `1e-8`
"""
__tolerance() = 1e-8

"""
$(TYPEDSIGNATURES)

Default maximum of iterations: `1000`
"""
__max_iterations() = 1000

"""
$(TYPEDSIGNATURES)

Default value for Ipopt print level: `5`
"""
__ipopt_print_level() = 5

"""
$(TYPEDSIGNATURES)

Default value for Ipopt mu strategy: `adaptive`
"""
__ipopt_mu_strategy() = "adaptive"

"""
$(TYPEDSIGNATURES)

Default value for Ipopt linear solver: `mumps`
"""
__ipopt_linear_solver() = "mumps"

#="""
$(TYPEDSIGNATURES)

Default value for MadNLP linear solver: `umfpack`
"""
__madnlp_linear_solver() = "umfpack"=#

"""
$(TYPEDSIGNATURES)

Default value for Knitro print level: `3`
"""
__knitro_print_level() = 3
