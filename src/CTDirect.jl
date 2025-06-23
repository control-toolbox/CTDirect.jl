module CTDirect

import CTModels

using DocStringExtensions
using SparseArrays

# includes
include("utils.jl")
include("default.jl")
include("docp.jl")
include("disc/common.jl")
include("disc/euler.jl")
include("disc/irk.jl")
include("disc/midpoint.jl")
include("disc/trapeze.jl")
include("solution.jl")
include("solve.jl")

end
