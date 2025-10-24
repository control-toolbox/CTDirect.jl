module CTDirectExtADNLP

using CTDirect
using DocStringExtensions
using ADNLPModels
using CTModels

"""
$(TYPEDSIGNATURES)

Build the NLP model for the DOCP (ADNLPModels version)

# Keyword arguments (optional)
* `show_time`: (:true, [:false]) show timing details from ADNLPModels
* `adnlp_backend`: backend for ADNLPModels ([`:optimized`], `:manual`, `:default`)
"""
function CTDirect.build_nlp!(
    docp::CTDirect.DOCP{
        <:CTDirect.Discretization,
        <:CTModels.Model,
        <:CTDirect.ADNLPBackend,
    },
    x0;
    adnlp_backend=CTDirect.__adnlp_backend(),
    show_time=false, #+default
    kwargs...,
)

    # redeclare objective and constraints functions
    f = x -> CTDirect.DOCP_objective(x, docp)
    c! = (c, x) -> CTDirect.DOCP_constraints!(c, x, docp)

    # unused backends (option excluded_backend = [:jprod_backend, :jtprod_backend, :hprod_backend, :ghjvprod_backend] does not seem to work)
    unused_backends = (
        hprod_backend=ADNLPModels.EmptyADbackend,
        jtprod_backend=ADNLPModels.EmptyADbackend,
        jprod_backend=ADNLPModels.EmptyADbackend,
        ghjvprod_backend=ADNLPModels.EmptyADbackend,
    )

    # call NLP problem constructor
    if adnlp_backend == :manual

        # build sparsity patterns for Jacobian and Hessian
        J_backend = ADNLPModels.SparseADJacobian(
            docp.dim_NLP_variables,
            f,
            docp.dim_NLP_constraints,
            c!,
            CTDirect.DOCP_Jacobian_pattern(docp),
        )
        H_backend = ADNLPModels.SparseReverseADHessian(
            docp.dim_NLP_variables,
            f,
            docp.dim_NLP_constraints,
            c!,
            CTDirect.DOCP_Hessian_pattern(docp),
        )
        backend_options = (
            gradient_backend=ADNLPModels.ReverseDiffADGradient,
            jacobian_backend=J_backend,
            hessian_backend=H_backend,
        )

    else
        # use backend preset
        backend_options = (backend=adnlp_backend,)
    end

    # build NLP
    nlp = ADNLPModel!(
        f,
        x0,
        docp.bounds.var_l,
        docp.bounds.var_u,
        c!,
        docp.bounds.con_l,
        docp.bounds.con_u;
        minimize=(!docp.flags.max),
        backend_options...,
        unused_backends...,
        show_time=show_time,
    )

    # set NLP in DOCP
    docp.nlp = nlp

    return nothing
end

end
