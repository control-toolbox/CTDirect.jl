# Beam example from bocop

@def beam begin
    t ∈ [ 0, 1 ], time
    x ∈ R², state
    u ∈ R, control
    x₁(0) == 0
    x₂(0) == 1
    x₁(1) == 0
    x₂(1) == -1
    ẋ(t) == [x₂(t), u(t)]
    0 ≤ x₁(t) ≤ 0.1    
    -10 ≤ u(t) ≤ 10
    ∫(u(t)^2) → min
end

return ((ocp=beam, obj=8.898598, name="beam"))
