# Build functional OCP solution from discrete DOCP solution

#+++ todo: add minimalist solution for new model
# with at least objective for testing purposes
"""
$(TYPEDSIGNATURES)
   
Build OCP functional solution from DOCP discrete solution (given as a SolverCore.GenericExecutionStats)
"""
struct OptimalControlSolution
    objective::Float64
    
    function OptimalControlSolution(docp::DOCP, docp_solution)
        if docp.is_maximization
            objective = -docp_solution.objective
        else
            objective = docp_solution.objective
        end
        return new(objective)
    end
end