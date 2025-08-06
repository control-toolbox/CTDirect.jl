module CTDirectExtExa

using CTDirect
using CTModels: CTModels

using DocStringExtensions

using ExaModels

"""
$(TYPEDSIGNATURES)

Build the NLP model for the DOCP (ExaModels version)

# Keyword arguments (optional)
* `grid_size`: number of time steps for the discretized problem ([250])
* `disc_method`: discretization method ([`:trapeze`], `:euler`)
* `exa_backend`: backend for ExaModels ([`nothing`])
"""
function CTDirect.build_nlp!(
    docp::CTDirect.DOCP,
    nlp_model::CTDirect.ExaBackend,
    x0;
    grid_size=CTDirect.__grid_size(),
    disc_method=CTDirect.__disc_method(),
    exa_backend=CTDirect.__exa_backend(),
    kwargs...,
)
    # check if MadNLPGPU is loaded when using GPU backend
    if !isnothing(exa_backend) && !CTDirect.package_loaded("MadNLPGPU")
        error("Please load MadNLPGU for ExaModels with CUDA: julia> using MadNLPGPU")
    end

    # build nlp
    # (time_grid != __time_grid()) || throw("non uniform time grid not available for nlp_model = :exa") # todo: remove when implemented in CTParser
    # set initial guess (NB. do not broadcast, apparently fails on GPU arrays)
    # NB unused final control in examodel / euler, hence the different x0 sizes
    build_exa = CTModels.get_build_examodel(docp.ocp)

    ocp = docp.ocp
    n = CTModels.state_dimension(ocp)
    m = CTModels.control_dimension(ocp)
    q = CTModels.variable_dimension(ocp)
    state = hcat([x0[(1 + i * (n + m)):(1 + i * (n + m) + n - 1)] for i in 0:grid_size]...) # grid_size + 1 states
    control = hcat(
        [
            x0[(n + 1 + i * (n + m)):(n + 1 + i * (n + m) + m - 1)] for
            i in 0:(grid_size - 1)
        ]...,
    ) # grid_size controls...
    control = [control control[:, end]] # ... todo: pass indeed to grid_size only for euler(_b), trapeze and midpoint
    variable = x0[(end - q + 1):end]
    docp.nlp, docp.exa_getter = build_exa(;
        grid_size=grid_size,
        backend=exa_backend,
        scheme=disc_method,
        init=(variable, state, control),
    )
    # remark: docp.nlp.meta.x0[1:docp.dim_NLP_variables] = -vcat(state..., control..., variable) # also work, and supersedes previous init via ExaModels start (itself overridden by init in solve)
    return nothing
end

end
