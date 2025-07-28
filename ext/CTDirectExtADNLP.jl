module CTDirectExtADNLP

using CTDirect

using DocStringExtensions

using ADNLPModels

"""
$(TYPEDSIGNATURES)

Build the NLP model for the DOCP (ADNLPModels version)

# Keyword arguments (optional)
* `show_time`: (:true, [:false]) show timing details from ADNLPModels
* `adnlp_backend`: backend for ADNLPModels ([`:optimized`], `:manual`, `:default`)
"""
function CTDirect.build_nlp(
    nlp_model::CTDirect.ADNLPBackend,
    docp::CTDirect.DOCP,
    x0;
    adnlp_backend=CTDirect.__adnlp_backend(),
    show_time=false, #+default
    matrix_free=false, #+default
    nlp_solver=nothing,
    kwargs...,
)

    # redeclare objective and constraints functions
    f = x -> CTDirect.DOCP_objective(x, docp)
    c! = (c, x) -> CTDirect.DOCP_constraints!(c, x, docp)

    # call NLP problem constructor
    if adnlp_backend == :manual

        # build sparsity pattern
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

        # build NLP with given patterns; disable unused backends according to solver info
        if (
            nlp_solver isa CTDirect.IpoptBackend ||
            nlp_solver isa CTDirect.MadNLPBackend ||
            nlp_solver isa CTDirect.KnitroBackend
        )
            nlp = ADNLPModel!(
                f,
                x0,
                docp.bounds.var_l,
                docp.bounds.var_u,
                c!,
                docp.bounds.con_l,
                docp.bounds.con_u;
                gradient_backend=ADNLPModels.ReverseDiffADGradient,
                jacobian_backend=J_backend,
                hessian_backend=H_backend,
                hprod_backend=ADNLPModels.EmptyADbackend,
                jtprod_backend=ADNLPModels.EmptyADbackend,
                jprod_backend=ADNLPModels.EmptyADbackend,
                ghjvprod_backend=ADNLPModels.EmptyADbackend,
                show_time=show_time,
                #excluded_backend = [:jprod_backend, :jtprod_backend, :hprod_backend, :ghjvprod_backend]
            )
        else
            nlp = ADNLPModel!(
                f,
                x0,
                docp.bounds.var_l,
                docp.bounds.var_u,
                c!,
                docp.bounds.con_l,
                docp.bounds.con_u;
                gradient_backend=ADNLPModels.ReverseDiffADGradient,
                jacobian_backend=J_backend,
                hessian_backend=H_backend,
                show_time=show_time,
            )
        end
    else
        # build NLP; disable unused backends according to solver info
        if (
            nlp_solver isa CTDirect.IpoptBackend ||
            nlp_solver isa CTDirect.MadNLPBackend ||
            nlp_solver isa CTDirect.KnitroBackend
        )
            nlp = ADNLPModel!(
                f,
                x0,
                docp.bounds.var_l,
                docp.bounds.var_u,
                c!,
                docp.bounds.con_l,
                docp.bounds.con_u;
                backend=adnlp_backend,
                hprod_backend=ADNLPModels.EmptyADbackend,
                jtprod_backend=ADNLPModels.EmptyADbackend,
                jprod_backend=ADNLPModels.EmptyADbackend,
                ghjvprod_backend=ADNLPModels.EmptyADbackend,
                show_time=show_time,
            )
        else
            # use manual settings including matrix_free
            nlp = ADNLPModel!(
                f,
                x0,
                docp.bounds.var_l,
                docp.bounds.var_u,
                c!,
                docp.bounds.con_l,
                docp.bounds.con_u;
                backend=adnlp_backend,
                show_time=show_time,
                matrix_free=matrix_free,
            )
        end
    end

    return nlp
end

end
