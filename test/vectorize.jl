function dynamics_vs(t, x, u)
    return [x[2], u]
end
function dynamics_ss(t, x, u)
    return x+u
end
function dynamics_sv(t, x, u)
    return x+u[1]-u[2]
end
function dynamics_vv(t, x, u)
    return [x[2], u[2]]
end

# remove additional x component too
function vectorize(fun, dimx, dimu)
    if dimx == 1
        if dimu == 1
            fun2 = (t, x, u) -> fun(t, x[1], u[1])
        else
            fun2 = (t, x, u) -> fun(t, x[1], u)
        end
    else
        if dimu == 1
            fun2 = (t, x, u) -> fun(t, x[1:dimx], u[1])
        else
            fun2 = (t, x, u) -> fun(t, x[1:dimx], u)
        end
    end
    return fun2
end

t = 0.
x = [1.,2.]
u = [0.,1.]

dynamics = vectorize(dynamics_ss, 1, 1)
println(dynamics(t, x, u))
dynamics = vectorize(dynamics_vs, 2, 1)
println(dynamics(t, x, u))
dynamics = vectorize(dynamics_sv, 1, 2)
println(dynamics(t, x, u))
dynamics = vectorize(dynamics_vv, 2, 2)
println(dynamics(t, x, u))
