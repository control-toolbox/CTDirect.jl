module CTDirectExt

using CTDirect
using CTBase
using DocStringExtensions

using JLD2
using JSON3

"""
$(TYPEDSIGNATURES)
  
Export OCP solution in JLD / JSON format
"""
function CTDirect.export_ocp_solution(sol::OptimalControlSolution; filename_prefix="solution", format=:JLD)
    if format == :JLD
        save_object(filename_prefix * ".jld2", sol)
    elseif format == :JSON
        blob = Dict(
            "objective" => sol.objective,
            "time_grid" => sol.time_grid,
            "state" => state_discretized(sol),
            "control" => control_discretized(sol),
            "costate" => costate_discretized(sol)[1:(end-1), :],
            "variable" => sol.variable,
        )
        open(filename_prefix * ".json", "w") do io
            JSON3.pretty(io, blob)
        end
    else
        error("Export_ocp_solution: unknow format (should be :JLD or :JSON): ", format)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
  
Read OCP solution in JLD / JSON format
"""
function CTDirect.import_ocp_solution(ocp::OptimalControlModel; filename_prefix="solution", format=:JLD)

    if format == :JLD
        return load_object(filename_prefix * ".jld2")
    elseif format == :JSON
        json_string = read(filename_prefix * ".json", String)
        blob = JSON3.read(json_string)

        # NB. convert vect{vect} to matrix
        return OptimalControlSolution(
            ocp,
            blob.time_grid,
            stack(blob.state, dims=1),
            stack(blob.control, dims=1),
            blob.variable,
            stack(blob.costate, dims=1);
            objective=blob.objective,
        )
    else
        error("Export_ocp_solution: unknow format (should be :JLD or :JSON): ", format)
    end
end


end
