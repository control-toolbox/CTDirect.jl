module CTDirectExt

    using CTDirect

    # load save export
    using JLD2
    
    """
    $(TYPEDSIGNATURES)
     
    Save OCP solution in JLD2 format
    """
    function CTDirect.save_OCP_solution(sol::OptimalControlSolution; filename_prefix="solution")
        save_object(filename_prefix * ".jld2", sol)
        return nothing
    end
        
    """
    $(TYPEDSIGNATURES)
     
    Load OCP solution in JLD2 format
    """
    function CTDirect.load_OCP_solution(filename_prefix="solution")
        return load_object(filename_prefix * ".jld2")
    end


end