println("testing: shape contract (1-D = scalar)")

# Load the scalar-style non-regression problem
if !isdefined(Main, :scalar_style_integrator)
    include("../problems/scalar_style_integrator.jl")
end

# --- type-recording fakes --------------------------------------------------
# Compute with `[1]`-indexing (guarded for the empty-control case), so the same fakes
# are valid whether the caller hands a Number or a length-1 vector -- only the
# *assertions* below distinguish before/after the coercion.
const SHAPE_SEEN = Dict{Symbol,Any}()

shape_rec_dyn!(r, t, x, u, v) = (
    SHAPE_SEEN[:dyn] = (x, u, v);
    r[1] = -x[1] + (isempty(u) ? 0.0 : u[1]);
    nothing
)
shape_rec_lag(t, x, u, v) = (
    SHAPE_SEEN[:lag] = (x, u, v);
    isempty(u) ? 0.0 : 0.5 * u[1]^2
)
shape_rec_may(x0, xf, v) = (SHAPE_SEEN[:may] = (x0, xf, v); 0.0)
shape_rec_path!(r, t, x, u, v) = (
    SHAPE_SEEN[:path] = (x, u, v);
    r[1] = x[1] + (isempty(u) ? 0.0 : u[1]);
    nothing
)
shape_rec_bnd!(r, x0, xf, v) = (
    SHAPE_SEEN[:bnd] = (x0, xf, v); r[1] = x0[1]; r[2] = xf[1]; nothing
)

const SHAPE_SCHEMES = (
    :euler, :euler_implicit, :trapeze, :midpoint,
    :gauss_legendre_2_constant_control, :gauss_legendre_2, :gauss_legendre_3,
)

# Build a recording OCP with declared dimensions n (state), m (control), q (variable).
function shape_build_recording(n::Int, m::Int, q::Int)
    pre = CTModels.PreModel()
    CTModels.time!(pre; t0=0.0, tf=1.0)
    CTModels.state!(pre, n)
    m > 0 && CTModels.control!(pre, m)
    q > 0 && CTModels.variable!(pre, q)

    CTModels.dynamics!(pre, shape_rec_dyn!)
    CTModels.objective!(pre, :min; lagrange=shape_rec_lag, mayer=shape_rec_may)
    CTModels.time_dependence!(pre; autonomous=false)

    CTModels.constraint!(pre, :path; f=shape_rec_path!, lb=[-Inf], ub=[Inf], label=:rec_path)
    CTModels.constraint!(
        pre, :boundary; f=shape_rec_bnd!, lb=[-Inf, -Inf], ub=[Inf, Inf], label=:rec_boundary
    )

    return CTModels.build(pre)
end

# Drive the recording OCP through __objective / __constraints! directly (no solver):
# enough to hit every user-function call boundary, and keeps this in the unit-test band.
function shape_drive(ocp, scheme::Symbol)
    empty!(SHAPE_SEEN)
    docp = CTDirect.DOCP(ocp, 3, 1, scheme, nothing)
    xu = 0.1 * ones(docp.dim_NLP_variables)
    CTDirect.__objective(xu, docp)
    CTDirect.__constraints!(zeros(docp.dim_NLP_constraints), xu, docp)
    return copy(SHAPE_SEEN)
end

@testset verbose = true showtiming = true "Shape contract: 1-D = scalar" begin

    for scheme in SHAPE_SCHEMES
        # --- CONTRACT: 1-D quantities reach the user functions as scalars ---
        @testset "$scheme (n=m=q=1)" begin
            seen = shape_drive(shape_build_recording(1, 1, 1), scheme)
            for k in (:dyn, :lag, :path)
                x, u, v = seen[k]
                @test x isa Number    # a ForwardDiff.Dual is still a Number
                @test u isa Number
                @test v isa Number
            end
            for k in (:may, :bnd)
                x0, xf, v = seen[k]
                @test x0 isa Number && xf isa Number && v isa Number
            end
        end

        # --- CONTRACT: n-D quantities stay vectors ---------------------------
        @testset "$scheme (n=m=q=2)" begin
            seen = shape_drive(shape_build_recording(2, 2, 2), scheme)
            x, u, v = seen[:dyn]
            @test x isa AbstractVector && length(x) == 2
            @test u isa AbstractVector && length(u) == 2
            @test v isa AbstractVector && length(v) == 2
        end

        # --- CONTRACT: mixed dimensions coerce independently -----------------
        @testset "$scheme (n=2, m=1, q=1)" begin
            x, u, v = shape_drive(shape_build_recording(2, 1, 1), scheme)[:dyn]
            @test x isa AbstractVector && length(x) == 2
            @test u isa Number
            @test v isa Number
        end
    end

    # --- CONTRACT: dim 0 is left alone (control-free problems) ---------------
    @testset "zero control dimension" begin
        u = shape_drive(shape_build_recording(2, 0, 0), :midpoint)[:dyn][2]
        @test u isa AbstractVector && isempty(u)
    end

    # --- UNIT: the coercion primitive -----------------------------------------
    @testset "Unit: _dim_coerce" begin
        @test CTDirect._dim_coerce(1) === only
        @test CTDirect._dim_coerce(2) === identity
        @test CTDirect._dim_coerce(0) === identity   # `only` would throw on the empty view
        @test only(2.0) === 2.0                      # idempotent on a scalar
        @test only([2.0]) === 2.0
    end

    # --- INTEGRATION: scalar-written user functions are now legal -------------
    @testset "scalar-style dynamics solves" begin
        prob = scalar_style_integrator()
        sol = solve_problem(prob; scheme=:midpoint, grid_size=50, display=false)
        @test CTModels.successful(sol)
        @test sol.objective >= 0

        # --- INTEGRATION: solution side honours 1-D = scalar --------------------
        @test CTModels.state(sol)(0.5) isa Number
        @test CTModels.control(sol)(0.5) isa Number
    end
end
