module CTDirect

import CTBase
import CTModels

using DocStringExtensions
using ADNLPModels               # docp model with AD
using LinearAlgebra             # norm and misc
using HSL
using SparseArrays

# includes
include("default.jl")
include("docp.jl")
include("disc/common.jl")
#include("disc/euler.jl")
include("disc/irk.jl")
include("disc/midpoint.jl")
include("disc/trapeze.jl")
include("solution.jl")
include("solve.jl")

# exports
export available_methods
export is_solvable
export direct_transcription
export set_initial_guess
export direct_solve
export build_OCP_solution

end
