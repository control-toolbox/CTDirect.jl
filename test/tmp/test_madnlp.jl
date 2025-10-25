using Pkg
Pkg.activate(joinpath(@__DIR__))

using NLPModels
using ExaModels
using MadNLPMumps
using ADNLPModels
using Ipopt
using NLPModelsIpopt
using DataFrames

# min f(x)
f(x) = 1+x[1]^2

function build_problem(problem_type::Symbol, sense::Symbol)
    minimize = sense == :min
    if problem_type == :exa
        c = ExaCore(; minimize=minimize)
        x = variable(c, 1; start=2.0)
        expr = minimize ? f(x) : -f(x)
        objective(c, expr)
        return ExaModel(c)
    elseif problem_type == :adnlp
        objective_function = minimize ? f : (x -> -f(x))
        return ADNLPModel(objective_function, [2.0]; minimize=minimize)
    else
        error("Unsupported problem type: $(problem_type)")
    end
end

function solve_problem(problem, solver_type::Symbol)
    if solver_type == :madnlp
        solver = MadNLPSolver(problem; print_level=MadNLP.ERROR)
        results = solve!(solver)
        println("typeof(results): $(typeof(results))")
        return results.objective
    elseif solver_type == :ipopt
        results = ipopt(problem; print_level=0)
        println("typeof(results): $(typeof(results))")
        return results.objective
    else
        error("Unsupported solver type: $(solver_type)")
    end
end

senses = [:min, :max]
solver_types = [:madnlp, :ipopt]
problem_types = [:exa, :adnlp]

results = DataFrame(;
    sense=String[], solver=String[], problem_type=String[], objective_value=Float64[]
)

for sense in senses
    for solver_type in solver_types
        for problem_type in problem_types
            problem = build_problem(problem_type, sense)
            objective_value = solve_problem(problem, solver_type)
            push!(
                results,
                (
                    sense=String(sense),
                    solver=String(solver_type),
                    problem_type=String(problem_type),
                    objective_value=objective_value,
                ),
            )
        end
    end
end

results
