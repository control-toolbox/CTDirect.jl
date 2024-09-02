# Olivier, tu as ici l'extension qui fait les sauvergardes / lectures de solution en format julia (jld2) et texte (json). Le load/save est assez trivial, mais on pourrait unifier avec le json en ajoutant un argument format=:jld [:json]. Cette extension irqit ensuite dans CTBase typiquement vu que l'aspect direct n'intervient pas

module CTDirectExt

using CTDirect
using CTBase
using DocStringExtensions

using JLD2  # load / save
using JSON3 # read / export

"""
$(TYPEDSIGNATURES)
  
Save OCP solution in JLD2 format
"""
function JLD2.save(sol::OptimalControlSolution; filename_prefix = "solution")
    save_object(filename_prefix * ".jld2", sol)
    return nothing
end

"""
$(TYPEDSIGNATURES)
     
Load OCP solution in JLD2 format
"""
function JLD2.load(filename_prefix = "solution")
    return load_object(filename_prefix * ".jld2")
end

"""
$(TYPEDSIGNATURES)
  
Export OCP solution in JSON format
"""
function CTDirect.export_ocp_solution(sol::OptimalControlSolution; filename_prefix = "solution")
    # fuse into save ?
    blob = Dict(
        "objective" => sol.objective,
        "time_grid" => sol.time_grid,
        "state" => state_discretized(sol),
        "control" => control_discretized(sol),
        "costate" => costate_discretized(sol)[1:(end - 1), :],
        "variable" => sol.variable,
    )
    open(filename_prefix * ".json", "w") do io
        JSON3.pretty(io, blob)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
  
Read OCP solution in JSON format
"""
function CTDirect.import_ocp_solution(ocp::OptimalControlModel; filename_prefix = "solution")
    # fuse into load ? 
    json_string = read(filename_prefix * ".json", String)
    blob = JSON3.read(json_string)

    # NB. convert vect{vect} to matrix
    return OptimalControlSolution(
        ocp,
        blob.time_grid,
        stack(blob.state, dims = 1),
        stack(blob.control, dims = 1),
        blob.variable,
        stack(blob.costate, dims = 1);
        objective = blob.objective,
    )
end

end
