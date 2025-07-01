module CTDirectExtExa

using CTDirect
import CTModels

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
function CTDirect.build_nlp(
    nlp_model::CTDirect.ExaBackend,
    docp::CTDirect.DOCP,
    x0;
    grid_size=CTDirect.__grid_size(),
    disc_method=CTDirect.__disc_method(),
    exa_backend=CTDirect.__exa_backend(),
    kwargs...,
)

    # check if MadNLPGPU is loaded when using GPU backend
    if exa_backend isa CUDABackend && !CTDirect.package_loaded("MadNLPGPU")
        error("Please load MadNLPGU for ExaModels with CUDA: julia> using MadNLPGPU")
    end

    # build nlp
    # debug: (time_grid != __time_grid()) || throw("non uniform time grid not available for nlp_model = :exa") # todo: remove when implemented in CTParser
    build_exa = CTModels.get_build_examodel(docp.ocp)
    nlp = build_exa(; grid_size = grid_size, backend = exa_backend, scheme = disc_method) 
    
    # set initial guess (NB. do not broadcast, apparently fails on GPU arrays)
    # NB unused final control in examodel / euler, hence the different x0 sizes
    nlp.meta.x0[1:docp.dim_NLP_variables] = x0  

    return nlp
end


end
