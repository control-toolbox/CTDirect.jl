################################################################
################################################################
# A integrer dans CTBase vu que c'est la sauvegarde d'une solution OCP
# copie de CTDirectExt.jl qui definit l'extension avec les weakdeps JLD2 / JSON3

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
function CTDirect.export_ocp_solution(sol::OptimalControlSolution; filename_prefix = "solution", format = :JLD)
    if format == :JLD
        save_object(filename_prefix * ".jld2", sol)
    elseif format == :JSON
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
    else
        error("Export_ocp_solution: unknow format (should be :JLD or :JSON): ", format)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
  
Read OCP solution in JLD / JSON format
"""
function CTDirect.import_ocp_solution(ocp::OptimalControlModel; filename_prefix = "solution", format = :JLD)

    if format == :JLD
        return load_object(filename_prefix * ".jld2")
    elseif format == :JSON 
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
    else
        error("Export_ocp_solution: unknow format (should be :JLD or :JSON): ", format)
    end        
end

end



################################################################
################################################################
# cf solution.jl: a mettre dans CTBase ? (aqua indique un type piracy car on ne passe pas d'argument avec des types specifiques a DOCP). Ce constructeur pourrait etre utilise pour toute solution sous forme discrete, pas seulement les methodes directes.
"""
$(TYPEDSIGNATURES)
    
Build OCP functional solution from DOCP vector solution (given as raw variables and multipliers plus some optional infos)
"""
function CTBase.OptimalControlSolution(
    ocp::OptimalControlModel,
    T,
    X,
    U,
    v,
    P;
    objective = 0,
    iterations = 0,
    constraints_violation = 0,
    message = "No msg",
    stopping = nothing,
    success = nothing,
    constraints_types = (nothing, nothing, nothing, nothing, nothing),
    constraints_mult = (nothing, nothing, nothing, nothing, nothing),
    box_multipliers = (nothing, nothing, nothing, nothing, nothing, nothing),
)
    dim_x = state_dimension(ocp)
    dim_u = control_dimension(ocp)
    dim_v = variable_dimension(ocp)

    # check that time grid is strictly increasing
    # if not proceed with list of indexes as time grid
    if !issorted(T, lt = <=)
        println(
            "WARNING: time grid at solution is not strictly increasing, replacing with list of indices...",
        )
        println(T)
        dim_NLP_steps = length(T) - 1
        T = LinRange(0, dim_NLP_steps, dim_NLP_steps + 1)
    end

    # variables: remove additional state for lagrange cost
    x = ctinterpolate(T, matrix2vec(X[:, 1:dim_x], 1))
    p = ctinterpolate(T[1:(end - 1)], matrix2vec(P[:, 1:dim_x], 1))
    u = ctinterpolate(T, matrix2vec(U[:, 1:dim_u], 1))

    # force scalar output when dimension is 1
    fx = (dim_x == 1) ? deepcopy(t -> x(t)[1]) : deepcopy(t -> x(t))
    fu = (dim_u == 1) ? deepcopy(t -> u(t)[1]) : deepcopy(t -> u(t))
    fp = (dim_x == 1) ? deepcopy(t -> p(t)[1]) : deepcopy(t -> p(t))
    var = (dim_v == 1) ? v[1] : v

    # misc infos
    infos = Dict{Symbol, Any}()
    infos[:constraints_violation] = constraints_violation

    # nonlinear constraints and multipliers
    control_constraints = t -> ctinterpolate(T, matrix2vec(constraints_types[1], 1))(t)
    mult_control_constraints = t -> ctinterpolate(T, matrix2vec(constraints_mult[1], 1))(t)
    state_constraints = t -> ctinterpolate(T, matrix2vec(constraints_types[2], 1))(t)
    mult_state_constraints = t -> ctinterpolate(T, matrix2vec(constraints_mult[2], 1))(t)
    mixed_constraints = t -> ctinterpolate(T, matrix2vec(constraints_types[3], 1))(t)
    mult_mixed_constraints = t -> ctinterpolate(T, matrix2vec(constraints_mult[3], 1))(t)

    # boundary and variable constraints
    boundary_constraints = constraints_types[4]
    mult_boundary_constraints = constraints_mult[4]
    variable_constraints = constraints_types[5]
    mult_variable_constraints = constraints_mult[5]

    # box constraints multipliers
    mult_state_box_lower = t -> ctinterpolate(T, matrix2vec(box_multipliers[1][:, 1:dim_x], 1))(t)
    mult_state_box_upper = t -> ctinterpolate(T, matrix2vec(box_multipliers[2][:, 1:dim_x], 1))
    mult_control_box_lower = t -> ctinterpolate(T, matrix2vec(box_multipliers[3][:, 1:dim_u], 1))(t)
    mult_control_box_upper = t -> ctinterpolate(T, matrix2vec(box_multipliers[4][:, 1:dim_u], 1))
    mult_variable_box_lower, mult_variable_box_upper = box_multipliers[5], box_multipliers[6]

    # build and return solution
    if is_variable_dependent(ocp)
        return OptimalControlSolution(
            ocp;
            state = fx,
            control = fu,
            objective = objective,
            costate = fp,
            time_grid = T,
            variable = var,
            iterations = iterations,
            stopping = stopping,
            message = message,
            success = success,
            infos = infos,
            control_constraints = control_constraints,
            state_constraints = state_constraints,
            mixed_constraints = mixed_constraints,
            boundary_constraints = boundary_constraints,
            variable_constraints = variable_constraints,
            mult_control_constraints = mult_control_constraints,
            mult_state_constraints = mult_state_constraints,
            mult_mixed_constraints = mult_mixed_constraints,
            mult_boundary_constraints = mult_boundary_constraints,
            mult_variable_constraints = mult_variable_constraints,
            mult_state_box_lower = mult_state_box_lower,
            mult_state_box_upper = mult_state_box_upper,
            mult_control_box_lower = mult_control_box_lower,
            mult_control_box_upper = mult_control_box_upper,
            mult_variable_box_lower = mult_variable_box_lower,
            mult_variable_box_upper = mult_variable_box_upper,
        )
    else
        return OptimalControlSolution(
            ocp;
            state = fx,
            control = fu,
            objective = objective,
            costate = fp,
            time_grid = T,
            iterations = iterations,
            stopping = stopping,
            message = message,
            success = success,
            infos = infos,
            control_constraints = control_constraints,
            state_constraints = state_constraints,
            mixed_constraints = mixed_constraints,
            boundary_constraints = boundary_constraints,
            mult_control_constraints = mult_control_constraints,
            mult_state_constraints = mult_state_constraints,
            mult_mixed_constraints = mult_mixed_constraints,
            mult_boundary_constraints = mult_boundary_constraints,
            mult_state_box_lower = mult_state_box_lower,
            mult_state_box_upper = mult_state_box_upper,
            mult_control_box_lower = mult_control_box_lower,
            mult_control_box_upper = mult_control_box_upper,
        )
    end
end
