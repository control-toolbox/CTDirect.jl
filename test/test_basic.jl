using CTDirect

# simple integrator min energy
ocp = Model()
state!(ocp, 1)
control!(ocp, 1)
time!(ocp, [0, 1])
constraint!(ocp, :initial, -1, :initial_constraint)
constraint!(ocp, :final, 0, :final_constraint)
dynamics!(ocp, (x, u) -> -x + u)
objective!(ocp, :lagrange, (x, u) -> u^2)

# all-in-one solve call
println("Test simple integrator: all in one solve call")
sol = solveDirect(ocp, grid_size=100, print_level=5, tol=1e-12)
println("Expected Objective 0.313, found ", sol.objective)

# split calls
println("Test simple integrator: split calls")
println("Direct transcription")
docp = directTranscription(ocp, grid_size=100)
nlp = getNLP(docp)
println("Solve discretized problem and retrieve solution")
sol = solveDOCP(docp, print_level=5, tol=1e-12)
println("Expected Objective 0.313, found ", sol.objective)

# fail test
#sol = solveDirect(ocp, grid_size=100, print_level=5, tol=1e-12, max_iter=1) ok

@def ocp2 begin
    tf ∈ R, variable
    t ∈ [ 0, tf ], time
    x ∈ R², state
    u ∈ R, control
    tf ≥ 0
    -1 ≤ u(t) ≤ 1
    q = x₁
    v = x₂
    q(0) == 1
    v(0) == 2
    q(tf) == 0
    v(tf) == 0
    0 ≤ q(t) ≤ 5
    -2 ≤ v(t) ≤ 3
    (u^2)(t) ≤ 100
    ẋ(t) == [ v(t), u(t) ]
    tf → min
end
sol = solveDirect(ocp2, grid_size=100, print_level=5, tol=1e-12)
