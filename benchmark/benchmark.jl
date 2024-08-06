# Benchmark
include("../test/deps.jl")
using Printf
import LinearAlgebra

#######################################################
# set environment

# choose NLP solver
nlp_solver = :madnlp

# linear solver: default mumps; spral, ma27, ma57, ma77, ma86, ma97
linear_solver = "mumps"
@printf("NLP Solver: %s   Linear solver: %s\n", nlp_solver, linear_solver)

# blas backend
using MKL # Replace OpenBLAS with Intel MKL
blas_config = LinearAlgebra.BLAS.lbt_get_config()
@printf("Blas config: %s\n", blas_config)
# AD backend ?


#######################################################
# set parameters
tol = 1e-8
grid_size = 1000
precompile = true
@printf("Settings: tol=%g grid_size=%d precompile=%s\n\n", tol, grid_size, precompile)

#######################################################
# load examples library
problem_path = pwd()*"/test/problems"
for problem_file in filter(contains(r".jl$"), readdir(problem_path; join=true))
    include(problem_file)
end

# load problems for benchmark
print("Loading problems: ")
names_list = ["beam", "bioreactor_1day_periodic", "fuller", "goddard", "jackson", "vanderpol"]
problem_list = []
for problem_name in names_list
    ocp_data = getfield(Main, Symbol(problem_name))()
    @printf("%s ", ocp_data.name)
    push!(problem_list,ocp_data)
end
println("")

#######################################################
# precompile if required
if precompile
    t_precomp = 0.
    print("Precompilation step: ")
    for problem in problem_list
        @printf("%s ",problem[:name])
        t = @elapsed solve(problem[:ocp], nlp_solver, linear_solver=linear_solver, max_iter=0, display=false)
        global t_precomp += t
    end
    @printf("\nPrecompilation total time %6.2f\n",t_precomp)
end

#######################################################
# solve examples with timer and objective check
t_list = []
println("Benchmark step")
for problem in problem_list
    t = @elapsed local sol = solve(problem[:ocp], nlp_solver, init=problem[:init], display=false, linear_solver=linear_solver, grid_size=grid_size, tol=tol)
    if !isnothing(problem[:obj]) && !isapprox(sol.objective, problem[:obj], rtol=5e-2)
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