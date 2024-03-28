# try to determine the efficiency of the scalar/vector handling

dimx = 3
dimu = 1
N = 10
xu = ones((N+1)*dimx+N*dimu)

function getx(xu, step)
    n = dimx
    if n == 1
        return xu[step*n+1]
    else
        return xu[step*n+1:step*n+n]
    end
end

function getu(xu, step)
    n = dimx
    m = dimu
    if m == 1
        return xu[(N+1)*n + step*m + 1]
    else
        return xu[(N+1)*n + step*m + 1 : (N+1)*n + step*m + m]
    end

function dynamics_scalar(x, u)
    f = zeros(dimx)
    f[1] = x[2]
    f[2] = u
    f[3] = u*u
    return f
end

function dynamics_vector(x, u)
    f = zeros(dimx)
    f[1] = x[2]
    f[2] = u[1]
    f[3] = u[1]*u[1]
    return f
end

#compile
step = 5
dynamics_scalar(getx(xu,step), getu(xu,step))
dynamics_vector(getx(xu,step), getu(xu,step))

#dynamic dispatch ?
step = 6
@code_warntype dynamics_scalar(getx(xu,step), getu(xu,step)) 
@code_warntype dynamics_vector(getx(xu,step), getu(xu,step)) 