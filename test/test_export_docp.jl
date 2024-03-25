using CTDirect

# simple integrator min energy
ocp1 = Model()
state!(ocp1, 1)
control!(ocp1, 1)
time!(ocp1, [0, 1])
constraint!(ocp1, :initial, -1, :initial_constraint)
constraint!(ocp1, :final, 0, :final_constraint)
dynamics!(ocp1, (x, u) -> -x + u)
objective!(ocp1, :lagrange, (x, u) -> u^2)

println("Test simple integrator: basic init")
docp = DirectTranscription(ocp1, grid_size=100)
nlp = getNLP(docp)
#print(nlp.meta.x0)
sol = solveDOCP(docp, print_level=5, tol=1e-12)
println("Expected Objective 0.313, found ", sol.objective)

# different starting guess
println("with constant init x=0.5 and u=0")
init_constant = OptimalControlInit(x_init=[-0.5], u_init=0)
setDOCPInit(docp, init_constant)
#print(docp.nlp.meta.x0)
sol = solveDOCP(docp, print_level=5, tol=1e-12)
println("Expected Objective 0.313, found ", sol.objective)

# init from solution
init_sol = OptimalControlInit(sol)
setDOCPInit(docp, init_sol)
#print(docp.nlp.meta.x0)
sol = solveDOCP(docp, print_level=5, tol=1e-12)
println("Expected Objective 0.313, found ", sol.objective)

# pass init directly to solve call
setDOCPInit(docp, OptimalControlInit()) # reset init in docp
sol = solveDOCP(docp, init=init_sol, print_level=5, tol=1e-12)
println("Expected Objective 0.313, found ", sol.objective)
#print(docp.nlp.meta.x0)
sol = solveDOCP(docp, print_level=5, tol=1e-12)
println("Expected Objective 0.313, found ", sol.objective)


# check types on objective and constraints functions
#@code_warntype ipopt_objective(xu, docp)
