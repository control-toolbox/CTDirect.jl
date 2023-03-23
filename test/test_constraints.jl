using CTDirect
using CTProblems
using CTBase # for the functions
using Test

println("direct_infos function tests")
# 
prob = Problem(:integrator, :dim2, :energy); 
ocp = prob.model
constraint!(ocp, :control, -4.01, 4.01, :control_con1)

init = [1., 0.5, 0.3]
nlp = CTDirect.ADNLProblem(ocp, 10, init)
println(nlp)
l_var, u_var = CTDirect.variables_bounds()
println("lvar = ", l_var)
println("uvar = ", u_var)

