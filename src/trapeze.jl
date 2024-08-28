

struct TrapezeTag <: DiscretizationTag 
    stage::Int
    additional_controls::Int 
    # add control at tf
    TrapezeTag() = new(0, 1)
end

