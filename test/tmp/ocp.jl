using Pkg
Pkg.activate(".")

# OptimalControl
using CTBase
using CTParser: CTParser, @def, prefix!, e_prefix!
using CTModels:
    CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTDirect:
    CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution, nlp
prefix!(:CTModels) # set CTParser def macro to use CTModels instead of OptimalControl
e_prefix!(:CTBase) # set CTParser def macro to use CTBase instead of OptimalControl (errors)

# activate NLP modelers
#using ADNLPModels
# + using ExaModels (in test_exa for now)

# activate NLP solvers
using NLPModelsIpopt
using MadNLP

@def ocp begin
    tf ∈ R, variable
    t ∈ [0, tf], time
    x ∈ R², state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    x(0) == [0, 0]
    x(tf) == [1, 0]
    0.05 ≤ tf ≤ Inf
    ẋ(t) == [x₂(t), u(t)]
    tf → min
end

sol = solve(ocp, :adnlp, :ipopt; print_level=5);
