module CTDirect

using ADNLPModels
using ExaModels
using CTBase
using CTModels: CTModels
using DocStringExtensions
using SparseArrays
using SolverCore: SolverCore
using NLPModels: NLPModels

# ----------------------------------------------------------------------
# TYPES

const AbstractOptimalControlProblem = CTModels.AbstractModel

# # ----------------------------------------------------------------------
# # EXTENSIONS

# # NLP solver backend extensions 
# abstract type AbstractNLPSolverBackend end
# struct IpoptBackend <: AbstractNLPSolverBackend end
# struct MadNLPBackend <: AbstractNLPSolverBackend end
# struct KnitroBackend <: AbstractNLPSolverBackend end

# # NLP model backend extensions
abstract type AbstractNLPModelBackend end
struct ADNLPBackend <: AbstractNLPModelBackend end
struct ExaBackend <: AbstractNLPModelBackend end

# ## Extensions and weak dependencies (see ext/CTDirectExt***)
# const WEAKDEPS = Dict{Type,Any}(
#     # NLP solver
#     IpoptBackend => [:NLPModelsIpopt],
#     MadNLPBackend => [:MadNLP],
#     KnitroBackend => [:NLPModelsKnitro],
#     # NLP modeller
#     ADNLPBackend => [:ADNLPModels],
#     ExaBackend => [:ExaModels],
# )

# includes
include("utils.jl")
include("default.jl")
include("core_types.jl")
include("discretization_api.jl")
include("collocation.jl")
include("collocation_core.jl")

include("disc/common.jl")
include("disc/euler.jl")
include("disc/irk.jl")
include("disc/midpoint.jl")
include("disc/trapeze.jl")

include("solution.jl")
include("solve.jl")

end
