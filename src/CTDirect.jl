module CTDirect

using CTBase
using DocStringExtensions
using NLPModelsIpopt, ADNLPModels # nlp modeling and resolution

# Other declarations
const nlp_constraints = CTBase.nlp_constraints
const __grid_size_direct() = 100
const __print_level_ipopt = CTBase.__print_level_ipopt
const __mu_strategy_ipopt = CTBase.__mu_strategy_ipopt
const __display = CTBase.__display
const matrix2vec = CTBase.matrix2vec

# includes
include("init.jl")
include("utils.jl")
include("problem.jl")
include("solution.jl")
include("solve.jl")

# export functions only for user
export solve
export is_solvable
export OptimalControlInit

end