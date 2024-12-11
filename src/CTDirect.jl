module CTDirect

using CTBase # later import CTBase and import CTModels
using DocStringExtensions
using ADNLPModels               # docp model with AD
using LinearAlgebra             # norm and misc
using HSL

import CTBase: OptimalControlSolution, CTBase   # extended

# other declarations
const matrix2vec = CTBase.matrix2vec

# includes
include("utils.jl")
include("default.jl")
include("docp.jl")
include("disc/midpoint.jl")
include("disc/trapeze.jl")
include("disc/irk.jl")
include("solution.jl")
include("solve.jl")

# exports
export available_methods
export is_solvable
export direct_transcription
export set_initial_guess
export direct_solve

end
