module CTDirect

using CTBase
using DocStringExtensions
using ADNLPModels               # docp model with AD
using NLPModelsIpopt            # NLP solver
using LinearAlgebra             # norm
using JLD2
using JSON3
using StructTypes

# Other declarations
const __grid_size_direct() = 100
const __print_level_ipopt = CTBase.__print_level_ipopt
const __mu_strategy_ipopt = CTBase.__mu_strategy_ipopt
const __display = CTBase.__display
const nlp_constraints = CTBase.nlp_constraints
const matrix2vec = CTBase.matrix2vec

# includes
include("utils.jl")
include("problem.jl")
include("solution.jl")
include("solve.jl")

# export functions only for user
export solve
export available_methods
export is_solvable
export directTranscription
export getNLP
export setDOCPInit
export OCPSolutionFromDOCP
export OCPSolutionFromDOCP_raw
export save_OCP_solution
export load_OCP_solution
export OCP_Solution_discrete


end