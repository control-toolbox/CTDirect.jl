# Benchmark
include("../test/deps.jl")
using Printf
import LinearAlgebra

#######################################################
# set environment
# linear solver: default mumps; spral, ma27, ma57, ma77, ma86, ma97
linear_solver = "mumps"
@printf("Profile: Linear solver: %s\n", linear_solver)
# blas backend
using MKL # Replace OpenBLAS with Intel MKL
blas_config = LinearAlgebra.BLAS.lbt_get_config()
@printf("Blas config: %s\n", blas_config)
# AD backend ?

#######################################################
# set parameters
tol = 1e-8
grid_size = 5000
precompile = true
@printf("Settings: tol=%g grid_size=%d precompile=%s\n\n", tol, grid_size, precompile)

#######################################################
# load examples
# +++ could add :init field in the problems !
print("Loading problems: ")
names_list = ["beam", "bioreactor_1day_periodic", "fuller", "goddard", "insurance", "jackson", "swimmer", "vanderpol"]
problem_list = []
problem_path = pwd()*"/problems"
for problem_file in filter(contains(r".jl$"), readdir(problem_path; join=true))
    ocp_data = include(problem_file)
    if ocp_data.name in names_list
        @printf("%s ", ocp_data.name)
        push!(problem_list,ocp_data)
    end
end
println("")

#######################################################
# precompile if required
if precompile
    t_precomp = 0.
    print("Precompilation step: ")
    for problem in problem_list
        @printf("%s ",problem[:name])
        t = @elapsed solve(problem[:ocp], linear_solver=linear_solver, max_iter=0, display=false)
        global t_precomp += t
    end
    @printf("\nPrecompilation total time %6.2f\n",t_precomp)
end

#######################################################
# solve examples with timer and objective check
t_list = []
println("Benchmark step")
for problem in problem_list
    t = @elapsed local sol = solve(problem[:ocp], display=false, linear_solver=linear_solver, grid_size=grid_size, tol=tol)
    if !isapprox(sol.objective, problem[:obj], rtol=5e-2)
        error("Objective mismatch for ", problem[:name], ": ", sol.objective, " instead of ", problem[:obj])
    else
        @printf("%-30s completed in %6.2f s after %4d iterations\n",problem[:name],t,sol.iterations)
        append!(t_list,t)
    end
end

#######################################################
# print total time
@printf("\nTotal time (s): %6.2f \n", sum(t_list))

# to redirect script to file (OK)
#=
using Suppressor
output = @capture_out include("benchmark/benchmark.jl")
open("log.txt", "w") do file
    write(file, output)
end
=#