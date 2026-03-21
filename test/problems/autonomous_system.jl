# Parameter estimation problems without control (dimension 0)
# Using CTModels directly since CTParser does not yet support zero control dimension
# CTModels must be loaded before including this file

# Problem 1: Estimate initial condition of harmonic oscillator
function estimate_initial_condition()
    pre = CTModels.PreModel()
    CTModels.time!(pre; t0=0.0, tf=π/2)  # quarter period
    CTModels.state!(pre, 2)
    CTModels.variable!(pre, 2)  # v = [x₁(0), x₂(0)] to estimate
    
    # Dynamics: harmonic oscillator ẋ = [-x₂, x₁]
    dynamics!(r, t, x, u, v) = r .= [-x[2], x[1]]
    CTModels.dynamics!(pre, dynamics!)
    
    # Objective: minimize distance to target final state [0, 1]
    # Solution should be v ≈ [1, 0]
    CTModels.objective!(pre, :min; 
        mayer=(x0, xf, v) -> (xf[1] - 0.0)^2 + (xf[2] - 1.0)^2
    )
    
    CTModels.definition!(pre, quote
        v ∈ R², variable
        t ∈ [0, π/2], time
        x ∈ R², state
        x(0) == v
        ẋ(t) == [-x₂(t), x₁(t)]
        (xf[1])^2 + (xf[2] - 1)^2 → min
    end)
    CTModels.time_dependence!(pre; autonomous=true)
    
    # Initial condition is the variable to estimate
    f_initial(r, x0, xf, v) = r .= x0 .- v
    CTModels.constraint!(pre, :boundary; 
        f=f_initial,
        lb=[0.0, 0.0], 
        ub=[0.0, 0.0],
        label=:initial
    )
    
    ocp = CTModels.build(pre)
    return ((ocp=ocp, obj=nothing, name="estimate_initial", init=()))
end

# Problem 2: Estimate parameter in dynamics (rotation rate)
function estimate_rotation_rate()
    pre = CTModels.PreModel()
    CTModels.time!(pre; t0=0.0, tf=1.0)
    CTModels.state!(pre, 2)
    CTModels.variable!(pre, 1)  # v = [α] rotation rate to estimate
    
    # Dynamics: ẋ = α*[-x₂, x₁] (rotation with unknown rate)
    dynamics!(r, t, x, u, v) = r .= v[1] .* [-x[2], x[1]]
    CTModels.dynamics!(pre, dynamics!)
    
    # Objective: minimize distance to target final state [0, 1]
    # and regularize parameter (solution should be α ≈ π/2)
    CTModels.objective!(pre, :min; 
        mayer=(x0, xf, v) -> (xf[1] - 0.0)^2 + (xf[2] - 1.0)^2 + 0.01*v[1]^2
    )
    
    CTModels.definition!(pre, quote
        α ∈ R, variable
        t ∈ [0, 1], time
        x ∈ R², state
        0 ≤ α ≤ 10
        x(0) == [1, 0]
        ẋ(t) == α * [-x₂(t), x₁(t)]
        (xf[1])^2 + (xf[2] - 1)^2 + 0.01*α^2 → min
    end)
    CTModels.time_dependence!(pre; autonomous=false)
    
    # Known initial condition
    f_initial(r, x0, xf, v) = r .= x0 .- [1.0, 0.0]
    CTModels.constraint!(pre, :boundary; 
        f=f_initial,
        lb=[0.0, 0.0], 
        ub=[0.0, 0.0],
        label=:initial
    )
    
    # Bounds on parameter
    CTModels.constraint!(pre, :variable; rg=1, lb=0.0, ub=10.0, label=:alpha_bounds)
    
    ocp = CTModels.build(pre)
    return ((ocp=ocp, obj=nothing, name="estimate_parameter", init=()))
end

# Problem 3: Least squares fit with path constraint
function least_squares_with_constraint()
    pre = CTModels.PreModel()
    CTModels.time!(pre; t0=0.0, tf=1.0)
    CTModels.state!(pre, 2)
    CTModels.variable!(pre, 2)  # v = [x₁(0), x₂(0)] to estimate
    
    # Dynamics: harmonic oscillator
    dynamics!(r, t, x, u, v) = r .= [-x[2], x[1]]
    CTModels.dynamics!(pre, dynamics!)
    
    # Objective: minimize distance to measurements along trajectory
    # Synthetic measurements at t=0.5: [0.7, 0.7]
    CTModels.objective!(pre, :min; 
        lagrange=(t, x, u, v) -> (t - 0.5)^2 * ((x[1] - 0.7)^2 + (x[2] - 0.7)^2),
        mayer=(x0, xf, v) -> 0.01 * (v[1]^2 + v[2]^2)
    )
    
    CTModels.definition!(pre, quote
        v ∈ R², variable
        t ∈ [0, 1], time
        x ∈ R², state
        x(0) == v
        ẋ(t) == [-x₂(t), x₁(t)]
        x₁(t)^2 + x₂(t)^2 ≤ 2
        ∫((t - 0.5)^2 * ((x₁(t) - 0.7)^2 + (x₂(t) - 0.7)^2)) + 0.01*(v[1]^2 + v[2]^2) → min
    end)
    CTModels.time_dependence!(pre; autonomous=true)
    
    # Initial condition is the variable
    f_initial(r, x0, xf, v) = r .= x0 .- v
    CTModels.constraint!(pre, :boundary; 
        f=f_initial,
        lb=[0.0, 0.0], 
        ub=[0.0, 0.0],
        label=:initial
    )
    
    # Path constraint: stay within radius
    f_path(r, t, x, u, v) = r .= x[1]^2 + x[2]^2
    CTModels.constraint!(pre, :path;
        f=f_path,
        lb=-Inf,
        ub=2.0,
        label=:radius
    )
    
    ocp = CTModels.build(pre)
    return ((ocp=ocp, obj=nothing, name="least_squares_constraint", init=()))
end
