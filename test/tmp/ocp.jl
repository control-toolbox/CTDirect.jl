using Pkg
Pkg.activate("dev-ctdirect"; shared=true)

# OptimalControl
using CTBase
using CTParser: CTParser, @def
using CTModels:
    CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTDirect:
    CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution

@def ocp begin
    tf ∈ R, variable
    t ∈ [0, tf], time
    x ∈ R², state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    x(0) == [0, 0]
    x(tf) == [1, 0]
    0.05 ≤ tf ≤ Inf
    ∂(x₁)(t) == x₂(t)
    ∂(x₂)(t) == u(t)
    tf → min
end

#
using ADNLPModels
using NLPModelsIpopt
sol = solve(ocp, :adnlp, :ipopt);

# 
using ExaModels
using MadNLP
#using MadNLPMumps
sol = solve(ocp, :exa, :madnlp);
