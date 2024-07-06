# Benchmark
include("common_deps.jl")
using Printf

# set environment
# default mumps
linear_solver = "spral"
@printf("Profile: linear_solver=%s\n\n", linear_solver)

# set parameters
tol = 1e-8
grid_size = 100
precompile = true
@printf("Settings: tol=%g grid_size=%d precompile=%s\n\n", tol, grid_size, precompile)

# load problems +++use loop later cf runtests
problem_list = []
push!(problem_list, include("problems/goddard.jl"))
push!(problem_list, include("problems/parametric.jl"))


# precompile if required
if precompile
    println("Precompilation step")
    for problem in problem_list
        @printf("%s ",problem[:name])
        solve(problem[:ocp], linear_solver=linear_solver, max_iter=1, display=false)
    end
    println("\n")
end

# solve problems with timer and objective check
t_list = []
println("Benchmark step")
for problem in problem_list
    t = @elapsed local sol = solve(problem[:ocp], display=false, linear_solver=linear_solver, grid_size=grid_size, tol=tol)
    if !isapprox(sol.objective, problem[:obj], rtol=1e-2)
        error("Objective mismatch for ", problem[:name], ": ", sol.objective, " instead of ", problem[:obj])
    else
        @printf("%-20s completed in %6.2f s\n",problem[:name],t)
        append!(t_list,t)
    end
end

# print total time
@printf("\nTotal time (s): %6.2f \n", sum(t_list))
