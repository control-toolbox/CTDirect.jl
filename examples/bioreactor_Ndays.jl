# Bioreactor example from bocop

# constants
beta = 1
c = 2
gamma = 1
halfperiod = 5
Ks = 0.05
mu2m = 0.1
mubar = 1
r = 0.005
pi = 3.141592653589793
T = 300

function growth(s)
    # MONOD
	return mu2m * s / (s+Ks)
end

function daynightgrowth(time)
    # light model: max^2 (0,sin) * mubar
	# DAY/NIGHT CYCLE: [0,2 halfperiod] rescaled to [0,2pi]
	days = time / (halfperiod*2)
    tau = (days - floor(days)) * 2*pi
    mu = max(0,sin(tau))^2 * mubar
	return mu
end

#= METHANE PROBLEM
mu2 according to growth model, mu according to light model
time scale is [0,10] for 24h (day then night)=#
@def bioreactor_N begin
    t ∈ [ 0, T ], time
    x ∈ R³, state
    u ∈ R, control
    y = x[1]
    s = x[2]
    b = x[3]
    mu = daynightgrowth(t)
    mu2 = growth(s(t))
    [0,0,0.001] ≤ x(t) ≤ [Inf, Inf, Inf]
    0 ≤ u(t) ≤ 1
    0.05 ≤ y(0) ≤ 0.25
    0.5 ≤ s(0) ≤ 5
    0.5 ≤ b(0) ≤ 3
    ẋ(t) == [mu*y(t)/(1+y(t))-(r+u(t))*y(t),
            -mu2*b(t) + u(t)*beta*(gamma*y(t)-s(t)),
            (mu2-u(t)*beta)*b(t)]
    ∫(mu2*b(t)/(beta+c)) → max
end

#sol=solve(bioreactor_N, grid_size=1000, init=(state=[50,50,50],))
#plot(sol)
return ((ocp=bioreactor_N, obj=19.0745, name="bioreactor_Ndays"))