# Goddard problem
# free final time, fixed final mass, max altitude

goddard = Model(variable=true)
Cd = 310
Tmax = 3.5
β = 500
b = 2
r0 = 1
v0 = 0
vmax = 0.1
m0 = 1
mf = 0.6
x0 = [ r0, v0, m0 ]
state!(goddard, 3)
control!(goddard, 1)
variable!(goddard, 1)
time!(goddard, t0=0, indf=1)
constraint!(goddard, :initial, lb=x0, ub=x0)
constraint!(goddard, :final, rg=3, lb=mf, ub=mf)
constraint!(goddard, :state, lb=[r0,v0,mf], ub=[r0+0.2,vmax,m0])
constraint!(goddard, :control, lb=0, ub=1)
constraint!(goddard, :variable, lb=0.01, ub=Inf)
objective!(goddard, :mayer,  (x0, xf, v) -> xf[1], :max)
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
dynamics!(goddard, (x, u, v) -> F0(x) + u*F1(x) )

# return problem and objective
return ((ocp=goddard, obj=1.01257, name="goddard"))