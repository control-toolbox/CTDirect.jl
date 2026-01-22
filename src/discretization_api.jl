function discretize(
    ocp::AbstractOptimalControlProblem, 
    discretizer::AbstractOptimalControlDiscretizer
)
    return discretizer(ocp)
end

__discretizer()::AbstractOptimalControlDiscretizer = Collocation()

function discretize(
    ocp::AbstractOptimalControlProblem;
    discretizer::AbstractOptimalControlDiscretizer=__discretizer(),
)
    return discretize(ocp, discretizer)
end