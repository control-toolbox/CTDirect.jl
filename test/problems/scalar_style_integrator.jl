# Scalar-style state/control problem: dynamics, lagrange cost and boundary function are
# written without any `[1]` indexing on their inputs. Illegal before the "1-D = scalar"
# call-boundary coercion (issue #610 §7 / #613), legal after.
# Built via CTModels' functional PreModel API rather than @def: the macro always emits
# `[1]`-indexed aliases regardless of the coercion, so it would not exercise this contract.
# CTModels must be loaded before including this file.

function scalar_style_integrator()
    pre = CTModels.PreModel()
    CTModels.time!(pre; t0=0.0, tf=1.0)
    CTModels.state!(pre, 1)
    CTModels.control!(pre, 1)

    # dynamics ẋ = -x + u -- x, u used directly, no [1]
    dynamics!(r, t, x, u, v) = (r[1] = -x + u; nothing)
    CTModels.dynamics!(pre, dynamics!)

    # running cost ∫ u^2 dt -- u used directly, no [1]
    CTModels.objective!(pre, :min; lagrange=(t, x, u, v) -> u^2)

    CTModels.time_dependence!(pre; autonomous=true)

    # boundary conditions x(0) = -1, x(1) = 0 -- x0, xf used directly, no [1]
    f_boundary(r, x0, xf, v) = (r[1] = x0 + 1.0; r[2] = xf; nothing)
    CTModels.constraint!(
        pre, :boundary; f=f_boundary, lb=[0.0, 0.0], ub=[0.0, 0.0], label=:endpoints
    )

    ocp = CTModels.build(pre)
    return ((ocp=ocp, obj=nothing, name="scalar_style_integrator", init=()))
end
