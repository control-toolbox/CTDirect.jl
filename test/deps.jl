using CTDirect
using CTBase

using NLPModelsIpopt
using MadNLP
using HSL

using JLD2
using JSON3

# local function for testing purposes only
# SHOULD BE A COPY OF THE ONE IN OPTIMALCONTROL ! 
# +++ later add a CTDirect.solve_direct intermediate function 
function solve(ocp::OptimalControlModel, description::Symbol...;
    init=nothing,
    grid_size::Int=CTDirect.__grid_size(),
    time_grid=CTDirect.__time_grid(),
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
    docp_solution = CTDirect.solve_docp(tag, docp, nlp; kwargs...)

    # build and return OCP solution
    return build_solution(docp, docp_solution)

 end
