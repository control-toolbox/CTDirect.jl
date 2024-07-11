# Fuller example

@def fuller begin
    t ∈ [ 0, 3.5 ], time
    x ∈ R², state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    x(0) == [ 0, 1 ]
    x(3.5) == [ 0, 0 ]
    ẋ(t) == [ x₂(t), u(t) ] 
    ∫(x₁(t)^2) → min
end

#sol = solve(fuller, grid_size=1000)
#plot(sol)

return((ocp=fuller, obj=2.683944e-1, name="fuller"))