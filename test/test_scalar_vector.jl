# try to determine the efficiency of the scalar/vector handling
n = 3
m = 1
N = 10
xu = ones((N+1)*n+N*m)

function getx(xu, step, n)
    if n == 1
        return xu[step*n+1]
    else
        return xu[step*n+1:step*n+n]
    end
end

function getu(xu, step, n, m)
    #if m == 1
    #    return xu[(N+1)*n + step*m + 1]
    #else
        return xu[(N+1)*n + step*m + 1 : (N+1)*n + step*m + m]
    #end
end

function dynamics_scalar(x, u, n)
    f = zeros(n)
    f[1] = x[2]
    f[2] = u
    f[3] = u*u
    return f
end

function dynamics_vector(x, u, n)
    f = zeros(n)
    f[1] = x[2]
    f[2] = u[1]
    f[3] = u[1]*u[1]
    return f
end

# ici u semble bien de type scalaire des la premiere evaluation
@code_warntype dynamics_scalar(getx(xu,5,3), getu(xu,5,3,1))
@code_warntype dynamics_vector(getx(xu,5,3), getu(xu,5,3,1))

# nb le code dynamics_scalar passe aussi avec un u vectoriel
# au niveau de CTBase qu'est ce qui coince exactement si on attend un argument scalaire et qu'on recupere un vecteur de taille 1 ?