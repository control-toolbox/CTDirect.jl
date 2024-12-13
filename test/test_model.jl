using CTDirect
import CTModels

# define problem with new model: simple integrator
pre_ocp = CTModels.PreModel()
CTModels.time!(pre_ocp, t0=0.0, tf=1.0)
CTModels.state!(pre_ocp, 1)
CTModels.control!(pre_ocp, 1)
CTModels.variable!(pre_ocp, 0)
f!(r, t, x, u, v) = r .= u .- x
CTModels.dynamics!(pre_ocp, f!)
l!(r, t, x, u, v) = r .= u .* u
CTModels.objective!(pre_ocp, :min, lagrange=l!)
bc!(r, x0, xf, v) = r .= [x0, xf]
CTModels.constraint!(pre_ocp, :boundary, f=bc!, lb=[-1, 0], ub=[-1, 0], label=:boundary)
CTModels.definition!(pre_ocp, Expr(:simple_integrator_min_energy))
ocp = CTModels.build_model(pre_ocp)

# direct transcription
docp,_ = direct_transcription(ocp)
