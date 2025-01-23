using CTDirect
using CTBase
using LinearAlgebra
using NLPModelsIpopt
using MKL
using BenchmarkTools
using Printf


#disc_method_list = [:gauss_legendre_2]
disc_method_list = [:trapeze, :gauss_legendre_2]
grid_size_list = [1000, 2000, 5000] 

# Jump
include("ab_jump.jl")
for disc_method in disc_method_list
    for grid_size in grid_size_list
        @printf("Jump %s %d:", disc_method, grid_size)
        @btime algal_bacterial_jump(grid_size=$grid_size, disc_method=$disc_method, print_level=0)
    end
end

# CTDirect
include("problems/algal_bacterial.jl")
for disc_method in disc_method_list
    for grid_size in grid_size_list
        @printf("CTDirect %s %d:", disc_method, grid_size)
        @btime direct_solve(algal_bacterial().ocp, grid_size=$grid_size, disc_method=$disc_method, print_level=0, constant_control=true)
    end
end
