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
function CTDirect.save(sol::OptimalControlSolution; filename_prefix="solution")
    save_object(filename_prefix * ".jld2", sol)
    return nothing
end
        
"""
$(TYPEDSIGNATURES)
     
Load OCP solution in JLD2 format
"""
function CTDirect.load(filename_prefix="solution")
    return load_object(filename_prefix * ".jld2")
end

"""
$(TYPEDSIGNATURES)
  
Export OCP solution in JSON format
"""
function CTDirect.export_OCP_solution(sol::OptimalControlSolution; filename_prefix="solution")
    open(filename_prefix * ".json", "w") do io
        JSON3.pretty(io, OCP_Solution_discrete(sol))
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
  
Read OCP solution in JSON format
"""
function CTDirect.import_OCP_solution(filename_prefix="solution")
    json_string = read(filename_prefix * ".json", String)
    return JSON3.read(json_string)
end


end