# Direct methods

"""
$(TYPEDSIGNATURES)

Used to set the default discretization method.
The default value is `trapeze`.
"""
__disc_method() = :trapeze

"""
$(TYPEDSIGNATURES)

Used to set the default grid size.
The default value is `250`.
"""
__grid_size() = 250

"""
$(TYPEDSIGNATURES)

Used to set the default initial guess.
The default value is `nothing` and will correspond to 0.1 for all variables.
"""
__ocp_init() = nothing

"""
$(TYPEDSIGNATURES)

Used to set the default display toggle.
The default value is true.
"""
__display() = true

"""
$(TYPEDSIGNATURES)

Used to set the default time grid.
The default value is `nothing`.
"""
__time_grid() = nothing


"""
$(TYPEDSIGNATURES)
Used to set the default backend for AD in ADNLPModels.
The default value is `:optimized`.
"""
__adnlp_backend() = :optimized

"""
$(TYPEDSIGNATURES)

Used to set the default tolerance.
The default value is `1e-8`.
"""
__tolerance() = 1e-8

"""
$(TYPEDSIGNATURES)

Used to set the default maximum of iterations.
The default value is `1000`.
"""
__max_iterations() = 1000

# IPOPT

"""
$(TYPEDSIGNATURES)

Used to set the default value of the print level of Ipopt for the direct method.
The default value is `5`.
"""
__ipopt_print_level() = 5

"""
$(TYPEDSIGNATURES)

Used to set the default value of the μ strategy of Ipopt for the direct method.
The default value is `adaptive`.
"""
__ipopt_mu_strategy() = "adaptive"

"""
$(TYPEDSIGNATURES)

Used to set the default value of the linear solver of Ipopt for the direct method.
The default value is `mumps`.
"""
__ipopt_linear_solver() = "mumps"

# MadNLP

"""
$(TYPEDSIGNATURES)

Used to set the default value of the linear solver of MadNLP for the direct method.
The default value is `umfpack`.
"""
__madnlp_linear_solver() = "umfpack"


"""
$(TYPEDSIGNATURES)

Used to set the default value of the print level of Knitro for the direct method.
The default value is `3`.
"""
__knitro_print_level() = 3
