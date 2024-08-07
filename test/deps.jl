using CTDirect
using CTBase

using NLPModelsIpopt
using MadNLP
using HSL

using JLD2
using JSON3

# local function for testing purposes only
# SHOULD BE A COPY OF THE ONE IN OPTIMALCONTROL ! 
# +++ some arguments are solver specific and should be defined only in the different calls to solve_docp ?
function solve(ocp::OptimalControlModel, description::Symbol...;
    init=CTDirect.__ocp_init_direct(),
    grid_size::Integer=CTDirect.__grid_size_direct(),
    display::Bool=CTDirect.__display(),
    print_level::Integer=CTDirect.__print_level_ipopt(),
    mu_strategy::String=CTDirect.__mu_strategy_ipopt(),
    max_iter::Integer=CTDirect.__max_iter(),
    tol::Real=CTDirect.__tol(),
    linear_solver::String=CTDirect.__linear_solver(),
    time_grid=nothing,
    kwargs...)

    method = getFullDescription(description, available_methods())
    #println(method)

    # build discretized OCP, including initial guess
    docp, nlp = direct_transcription(ocp, description, init=init, grid_size=grid_size, time_grid=time_grid)

    # solve DOCP
    if :ipopt ∈ method
        tag = CTDirect.IpoptTag()
    elseif :madnlp ∈ method
        tag = CTDirect.MadNLPTag()
    else
        error("no known solver in method", method)
    end
    docp_solution = CTDirect.solve_docp(tag, docp, nlp, display=display, tol=tol, max_iter=max_iter, kwargs...)

    # build and return OCP solution
    return build_solution(docp, docp_solution)

 end
