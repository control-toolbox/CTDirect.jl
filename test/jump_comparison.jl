using CTDirect: CTDirect, solve, direct_transcription, set_initial_guess, build_OCP_solution
using CTModels: CTModels, objective, state, control, variable, costate, time_grid, iterations
using CTParser: CTParser, @def, set_prefix
set_prefix(:CTModels) # tell CTParser def macro to use CTModels instead of OptimalControl

using LinearAlgebra
using NLPModelsIpopt
using MKL
using BenchmarkTools
using Printf

jump = false
ctdirect = true
#adnlp_backend_list = [:manual]
#adnlp_backend_list = [:optimized]
adnlp_backend_list = [:manual, :optimized]
#disc_method_list = [:gauss_legendre_2]
#disc_method_list = [:trapeze]
disc_method_list = [:trapeze, :gauss_legendre_2]
grid_size_list = [1000, 2000, 5000]

# Jump
include("jump/algal_bacterial_jump.jl")
if jump
    for disc_method in disc_method_list
        for grid_size in grid_size_list
            @printf("Jump %s %d:", disc_method, grid_size)
            @btime algal_bacterial_jump(grid_size=($grid_size), disc_method=($disc_method), print_level=0)
        end
    end
end

# CTDirect
include("problems/algal_bacterial.jl")
if ctdirect
    for backend in adnlp_backend_list
        for disc_method in disc_method_list
            for grid_size in grid_size_list
                @printf("CTDirect (%s) %s %d:", backend, disc_method, grid_size)
                @btime solve(algal_bacterial().ocp, tol=1e-8, grid_size=($grid_size), disc_method=($disc_method), print_level=0, adnlp_backend=($backend))
            end
        end
    end
end
