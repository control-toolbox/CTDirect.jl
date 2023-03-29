module CTDirect
using DocStringExtensions

# using
#
using CTBase

# nlp modeling and resolution
using NLPModelsIpopt, ADNLPModels

# Other declarations
const nlp_constraints = CTBase.nlp_constraints
const __grid_size_direct = CTBase.__grid_size_direct
const __print_level_ipopt = CTBase.__print_level_ipopt
const __mu_strategy_ipopt = CTBase.__mu_strategy_ipopt
const __display = CTBase.__display
const matrix2vec = CTBase.matrix2vec

# includes
include("utils.jl")
include("problem.jl")
include("solution.jl")
include("solve.jl")

# export functions only for user
export solve

end