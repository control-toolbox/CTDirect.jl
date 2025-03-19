module CTDirect

import CTBase
import CTModels

using ADNLPModels

using DocStringExtensions
using SparseArrays

# includes
include("default.jl")
include("docp.jl")
include("disc/common.jl")
include("disc/euler.jl")
include("disc/irk.jl")
include("disc/midpoint.jl")
include("disc/trapeze.jl")
include("solution.jl")
include("solve.jl")

# exports
export direct_transcription
export set_initial_guess
export solve

end
