module CTDirect

using CTBase
using DocStringExtensions
using ADNLPModels               # docp model with AD
using LinearAlgebra             # norm

# Other declarations
const nlp_constraints! = CTBase.nlp_constraints!
const matrix2vec = CTBase.matrix2vec

# includes
include("utils.jl")
include("default.jl")
include("problem.jl")
include("solution.jl")
include("solve.jl")

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