module CTDirect

using CTBase
using DocStringExtensions
using ADNLPModels               # docp model with AD
using LinearAlgebra             # norm
using CommonSolve: solve

# Other declarations (+++ sort among modules)
const __grid_size_direct() = 100
const __print_level_ipopt = CTBase.__print_level_ipopt
const __mu_strategy_ipopt = CTBase.__mu_strategy_ipopt
const __tol() = 1e-8
const __max_iter() = 1000
const __display = CTBase.__display
const nlp_constraints! = CTBase.nlp_constraints!
const matrix2vec = CTBase.matrix2vec
const __linear_solver() = "ma57"

# includes
include("utils.jl")
include("problem.jl")
include("solution.jl")
include("solve.jl")

# re exports
#export solve # from CommonSolve, extended in CTSolveExt

# exports
export available_methods
export is_solvable
export direct_transcription
export get_nlp
export set_initial_guess
export build_solution
export save
export load
export export_ocp_solution
export import_ocp_solution
export DOCP
export OCPDiscreteSolution

end