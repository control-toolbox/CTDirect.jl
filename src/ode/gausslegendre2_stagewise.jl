# ============================================================================
# src/ode/irk_gl2_stagewise.jl
#
# Minimal specialized stagewise version of Gauss-Legendre 2 collocation.
#
# Internal layout for NLP variables (vecteur xu):
# [X_0, U_0^1, U_0^2, K_0^1, K_0^2,
#  X_1, U_1^1, U_1^2, K_1^1, K_1^2,
#  ...
#  X_N-1, U_N-1^1, U_N-1^2, K_N-1^1, K_N-1^2,
#  X_N, V]
#
#
#
#
#
#
# ============================================================================

abstract type GenericIRKStagewise <: Scheme end

"""
Implicit Gauss-Legendre collocation for s=2 with stagewise controls.
"""
struct Gauss_Legendre_2_Stagewise <: GenericIRKStagewise
    info::String
    stage::Int
    butcher_a::Matrix{Float64}
    butcher_b::Vector{Float64}
    butcher_c::Vector{Float64}
    _control_block::Int
    _step_variables_block::Int
    _state_stage_eqs_block::Int
    _step_pathcons_block::Int
    _final_control::Bool

    function Gauss_Legendre_2_Stagewise(dims::DOCPdims, time::DOCPtime)
        control_block,
        step_variables_block,
        state_stage_eqs_block,
        step_pathcons_block,
        dim_NLP_variables,
        dim_NLP_constraints = GL2_stagewise_dims(dims, time)

        disc = new(
            "Implicit Gauss-Legendre collocation for s=2, 4th order, stagewise controls",
            2,
            [0.25 (0.25 - sqrt(3) / 6);
             (0.25 + sqrt(3) / 6) 0.25],
            [0.5, 0.5],
            [0.5 - sqrt(3) / 6, 0.5 + sqrt(3) / 6],
            control_block,
            step_variables_block,
            state_stage_eqs_block,
            step_pathcons_block,
            false,
        )

        return disc, dim_NLP_variables, dim_NLP_constraints
    end
end

"""
Dimensions for the specialized Gauss-Legendre-2 stagewise scheme.
    rappel : 
    dims.NLP_x : nombre de variables d’état
    dims.NLP_u : nombre de variables de contrôle
    time.control_steps : pas bien compris je suppose à 1 ici 
    dims.NLP_v : nombre de variables d’opti supplémentaires  (temps final ou initial libre par ex)

    retour : 
    step_variables_block : contient le nombre de variable à mémoriser pour un seul timestep 
    state_stage_eqs_block : combien size de vecteur d'état on stocke pour les équation ? equ d'état et les deux stages 
    dim_NLP_variables : le nombre total de variables d’optimisation du problème , tous les pas plus le pas final 
    dim_NLP_constraints : pour chaque pas de temps*((les contraintes d’état + stage)+ les contraintes de chemin) + les contraintes de chemin au temps final + les contraintes de bord

"""
function GL2_stagewise_dims(dims::DOCPdims, time::DOCPtime)
    stage = 2

    control_block = dims.NLP_u * time.control_steps

    # per step:
    # x_i + U_i^1 + U_i^2 + K_i^1 + K_i^2
    step_variables_block = dims.NLP_x + 2 * control_block + stage * dims.NLP_x

    # state equation + 2 stage equations
    state_stage_eqs_block = dims.NLP_x * (1 + stage)

    # same convention as irk.jl: path constraints at time steps
    step_pathcons_block = dims.path_cons

    dim_NLP_variables = time.steps * step_variables_block + dims.NLP_x + dims.NLP_v

    dim_NLP_constraints =
        time.steps * (state_stage_eqs_block + step_pathcons_block) +
        step_pathcons_block +
        dims.boundary_cons

    return control_block,
           step_variables_block,
           state_stage_eqs_block,
           step_pathcons_block,
           dim_NLP_variables,
           dim_NLP_constraints
end

# ----------------------------------------------------------------------------
# Indexing helpers
# ----------------------------------------------------------------------------

"""
Return control at stage j on time step i.
Convention: 1 <= i <= time.steps, j in {1,2}.
"""
function get_stagecontrol_at_time_step(xu, docp::DOCP{<:Gauss_Legendre_2_Stagewise}, i::Int, j::Int)
    disc = disc_model(docp)
    dims = docp.dims

    @assert 1 <= j <= 2
    #decalage pour se mettre au départ du block 
    #X_i, U_i^1, U_i^2, K_i^1, K_i^2,
    var_offset = (i - 1) * disc._step_variables_block
    #on se retrouve la :  U_i^1, U_i^2, K_i^1, K_i^2
    x_block_end = var_offset + dims.NLP_x
    #on déf départ et fin en fct de la taille du controle 
    ctrl0 = x_block_end + (j - 1) * disc._control_block + 1
    ctrl1 = x_block_end + j * disc._control_block
    #on retourne la view 
    return @view xu[ctrl0:ctrl1]
end

"""
Compatibility accessor: return a "control at time step".
By default we return U_i^1.
On retourne le controle N à N+1 si j n'est pas précisé on 
retourne le controle j=1 
"""
#gestion du dernier pas de temps 
function get_OCP_control_at_time_step(
    xu,
    docp::DOCP{<:Gauss_Legendre_2_Stagewise},
    i::Int,
    j::Int,
)
    if i == docp.time.steps + 1
        i = docp.time.steps
    end

    return get_stagecontrol_at_time_step(xu, docp, i, j)
end
#multiple dispatch si j n'est pas précisé 
function get_OCP_control_at_time_step(xu, docp::DOCP{<:Gauss_Legendre_2_Stagewise}, i::Int)
    return get_OCP_control_at_time_step(xu, docp, i, 1)
end

"""
Return stage variables K_i^j on time step i.
Convention: 1 <= i <= time.steps, j in {1,2}.
On extrait K_i^j
"""
function get_stagevars_at_time_step(xu, docp::DOCP{<:Gauss_Legendre_2_Stagewise}, i::Int, j::Int)
    disc = disc_model(docp)
    dims = docp.dims

    @assert 1 <= j <= 2

    var_offset = (i - 1) * disc._step_variables_block
    stage_block_start = var_offset + dims.NLP_x + 2 * disc._control_block
    k0 = stage_block_start + (j - 1) * dims.NLP_x + 1
    k1 = stage_block_start + j * dims.NLP_x

    return @view xu[k0:k1]
end

# ----------------------------------------------------------------------------
# Work array
# ----------------------------------------------------------------------------

"""
Set work array for dynamics and lagrange cost evaluations.
Layout:
- [x_i^j ; sum_bk]
"""
function setWorkArray(docp::DOCP{<:Gauss_Legendre_2_Stagewise}, xu, time_grid, v)
    dims = docp.dims
    work = similar(xu, dims.NLP_x + dims.NLP_x)
    return work
end

# ----------------------------------------------------------------------------
# Bounds + initial guess
# ----------------------------------------------------------------------------
"""
Donner les bornes de chacune des variables du NLP 
var_l <= xu <= var_u
var_l et var_u ont la meme structure que xu donc on peut 
réutiliser les meme fonction pour accéder au bon endroit
"""
function __variables_bounds!(docp::DOCP{<:Gauss_Legendre_2_Stagewise})
    #on récupère les bornes et le modele continu 
    var_l = docp.bounds.var_l
    var_u = docp.bounds.var_u
    ocp = ocp_model(docp)
    #
    x_lb, x_ub = build_bounds_block(
        docp.dims.NLP_x,
        CTModels.dim_state_constraints_box(ocp),
        CTModels.state_constraints_box(ocp),
    )
    u_lb, u_ub = build_bounds_block(
        docp.dims.NLP_u,
        CTModels.dim_control_constraints_box(ocp),
        CTModels.control_constraints_box(ocp),
    )
    #we change var_l and var_u at every time step with the 
    #right bound 
    for i in 1:(docp.time.steps + 1)
        set_state_at_time_step!(var_l, x_lb, docp, i)
        set_state_at_time_step!(var_u, x_ub, docp, i)
    end
    #si il y a une variable de controle (modif nouvelle version controle 0)
    if docp.dims.NLP_u > 0
        #for each timestep 
        for i in 1:docp.time.steps
            #for each control stage 
            for j in 1:2
                # on accède au bon indice et on utilise la view pour Set 
                # à modifier avec un setter propre 
                get_stagecontrol_at_time_step(var_l, docp, i, j) .= u_lb
                get_stagecontrol_at_time_step(var_u, docp, i, j) .= u_ub
            end
        end
    end
    #same thing for supplementary variables 
    if docp.dims.NLP_v > 0
        v_lb, v_ub = build_bounds_block(
            docp.dims.NLP_v,
            CTModels.dim_variable_constraints_box(ocp),
            CTModels.variable_constraints_box(ocp),
        )
        #stocked at the end of xu so at the end of var_l and var_u 
        #too so old setter are ok 
        set_optim_variable!(var_l, v_lb, docp)
        set_optim_variable!(var_u, v_ub, docp)
    end

    return var_l, var_u
end

function __initial_guess(
    docp::DOCP{<:Gauss_Legendre_2_Stagewise},
    init::CTModels.InitialGuess,
)   
    #everything to 0.1 
    NLP_X = 0.1 * ones(docp.dim_NLP_variables)
    disc = disc_model(docp)
    #si on a des init variable on mets les valeurs init dispo à la bonne place 
    set_optim_variable!(NLP_X, init.variable, docp)
    #on remplie les ti avec init state 
    time_grid = get_time_grid(NLP_X, docp)
    for i in 1:(docp.time.steps + 1)
        ti = time_grid[i]
        set_state_at_time_step!(NLP_X, init.state(ti), docp, i)
    end
    # si controle 
    if docp.dims.NLP_u > 0
        for i in 1:docp.time.steps
            ti = time_grid[i]
            hi = time_grid[i + 1] - ti

            for j in 1:disc.stage
                tij = ti + disc.butcher_c[j] * hi
                uij = init.control(tij)

                if !isnothing(uij)
                    get_stagecontrol_at_time_step(NLP_X, docp, i, j) .= uij
                end
            end
        end
    end

    return NLP_X
end

# ----------------------------------------------------------------------------
# Running cost
# ----------------------------------------------------------------------------

"""
Compute running cost using Gauss quadrature and stagewise controls.
"""
function integral(docp::DOCP{<:Gauss_Legendre_2_Stagewise}, xu, v, time_grid, f)
    disc = disc_model(docp)
    dims = docp.dims

    value = 0.0
    work_xij = similar(xu, dims.NLP_x)
    #pour chaque pas de temps i 
    for i in 1:docp.time.steps
        #temps de départ ti 
        ti = time_grid[i]
        #prendre l’état au début xi
        xi = get_OCP_state_at_time_step(xu, docp, i)
        #calculer la taille du pas hi
        hi = time_grid[i + 1] - ti

        local_sum = 0.0
        #pour chaque stage 
        for j in 1:2
            #calculer le temps du stage
            tij = ti + disc.butcher_c[j] * hi
            #récupérer le contrôle de stage uij
            uij = get_stagecontrol_at_time_step(xu, docp, i, j)
            #reconstruire l’état au stage xij
            #on copie x_i dans work_xij
            @views @. work_xij = xi
            #pour chaque stage (implicite)
            for l in 1:2
                #Récupération de K_i^l
                kil = get_stagevars_at_time_step(xu, docp, i, l)
                #Mise à jour de work_xij
                @views @. work_xij = work_xij + hi * disc.butcher_a[j, l] * kil
            end
            #evaluer le cout en f(tij,xij,uij,v) et l'ajouter pondéré par bj 
            local_sum += disc.butcher_b[j] * f(tij, work_xij, uij, v)
        end
        #ajouter hi * somme_locale à la valeur totale
        value += hi * local_sum
    end

    return value
end

# ----------------------------------------------------------------------------
# Dynamics + stage constraints
# ----------------------------------------------------------------------------

"""
Set state and stage constraints.
Convention: 1 <= i <= docp.time.steps
"""
function stepStateConstraints!(
    docp::DOCP{<:Gauss_Legendre_2_Stagewise},
    c,
    xu,
    v,
    time_grid,
    i,
    work,
)
    ocp = ocp_model(docp)
    disc = disc_model(docp)
    dims = docp.dims

    # work layout: [x_ij ; sum_bk]
    #allocate for work defined in setWorkArray
    work_xij = @view work[1:dims.NLP_x]
    work_sumbk = @view work[(dims.NLP_x + 1):(2 * dims.NLP_x)]
    #calculer l'endroit ou commencent les contraintes du pas i 
    offset = (i - 1) * (disc._state_stage_eqs_block + disc._step_pathcons_block)
    #début fin et taille du pas 
    ti = time_grid[i]
    tip1 = time_grid[i + 1]
    hi = tip1 - ti

    #etat aux extremités du pas 
    xi = get_OCP_state_at_time_step(xu, docp, i)
    xip1 = get_OCP_state_at_time_step(xu, docp, i + 1)
    #première dim réservées à équation d'état 
    offset_stage_eqs = dims.NLP_x
    #pour les stages 
    for j in 1:2
        tij = ti + disc.butcher_c[j] * hi
        kij = get_stagevars_at_time_step(xu, docp, i, j)
        uij = get_stagecontrol_at_time_step(xu, docp, i, j)
        #construire Σ b_j K_i^j
        if j == 1
            @views @. work_sumbk = disc.butcher_b[j] * kij
        else
            @views @. work_sumbk = work_sumbk + disc.butcher_b[j] * kij
        end
        #x_i^j = x_i + h_i Σ_l a_{j,l} K_i^l
        @views @. work_xij = xi
        for l in 1:2
            kil = get_stagevars_at_time_step(xu, docp, i, l)
            @views @. work_xij = work_xij + hi * disc.butcher_a[j, l] * kil
        end
        #on calcule f(t_i^j, x_i^j, u_i^j, v)
        CTModels.dynamics(ocp)(
            @view(c[(offset + offset_stage_eqs + 1):(offset + offset_stage_eqs + dims.NLP_x)]),
            tij,
            work_xij,
            uij,
            v,
        )
        #on stocke K_i^j - f(t_i^j, x_i^j, u_i^j, v)
        # rappel @views @. -> sans réallouer élément par élément 
        @views @. c[(offset + offset_stage_eqs + 1):(offset + offset_stage_eqs + dims.NLP_x)] =
            kij - c[(offset + offset_stage_eqs + 1):(offset + offset_stage_eqs + dims.NLP_x)]
        #offset contrainte suivante 
        offset_stage_eqs += dims.NLP_x
    end
    #contrainte d'évolution de l'état x_{i+1} - (x_i + h_i Σ_j b_j K_i^j)
    @views @. c[(offset + 1):(offset + dims.NLP_x)] =
        xip1 - (xi + hi * work_sumbk)

    return nothing
end
