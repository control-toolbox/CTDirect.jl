using CTDirect
using CTBase
using Test


n = 1
m = 1
t0 = 0
tf = 1
x0 = 2
xf = 0
ocp = Model()
state!(ocp, n)
control!(ocp, m)
time!(ocp, [t0, tf])
constraint!(ocp, :initial, x0, :initial_constraint)
constraint!(ocp, :control, -1, 1, :control_constraint)
dynamics!(ocp, (x, u) -> -x + u)
objective!(ocp, :lagrange, (x, u) -> x[1]^2 + u^2)
f = ocp.dynamics
lag1 = ocp.lagrange
v = Real[]

function test_scalar_vector(f)
    u_vector = false
    x_vector = false
    v = Real[]
    if n == 1 && m == 1
        x = zeros(n)
        u = zeros(m)
        try
            f(0, x, u, v)
            u_vector = true
            x_vector = true
        catch

        end
    end
    return x_vector, u_vector
end

x_vector, u_vector = test_scalar_vector(f)
println("u_vector = ", u_vector)
println("x_vector = ", x_vector)

@testset verbose = true showtiming = true "x and u scalar " begin
    @test f(t0, 2, 1, v) == -1          # pass
    @test f(t0, [1, 2], [1, 1], v) == [0, -1] # pass
    #   @test f(t0,2,[1],v) == -1        # fail
    #   @test f(t0,[2],1,v) == -1        # fail
    #   @test f(t0,[2],[1],v) == -1      # fail
    @test f(t0, [2], [1], v) == [-1]      # pass
    #   @test f(t0,[2,4],[1],v) == -1    # fail
    #   @test f(t0,x0[1:1],[1],v) == -1  # fail
    #   @test lag(t0,x0[1],1,v) == 5      # fail
end

n = 1
m = 1
t0 = 0
tf = 1
x0 = [2]
xf = [0]
ocp = Model()
state!(ocp, n)
control!(ocp, m)
time!(ocp, [t0, tf])
constraint!(ocp, :initial, x0, :initial_constraint)
constraint!(ocp, :control, -1, 1, :control_constraint)
dynamics!(ocp, (x, u) -> -x[1] + u[1])
objective!(ocp, :lagrange, (x, u) -> x[1]^2 + u^2)
f_vec = ocp.dynamics
v = Real[]
f_vec(t0, x0, 1, v)
lag_scal_vec = ocp.lagrange
@testset verbose = true showtiming = true "x and u vector of dim 1 " begin
    @test f_vec(t0, 2, 1, v) == -1
    @test f_vec(t0, 2, [1], v) == -1
    @test f_vec(t0, [2], 1, v) == -1
    @test f_vec(t0, [2], [1], v) == -1
    @test f_vec(t0, [2, 4], [1], v) == -1
    @test f_vec(t0, x0[1:1], [1], v) == -1
    @test lag_scal_vec(t0, x0[1:1], 1, v) == 5
end
