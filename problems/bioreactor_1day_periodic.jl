# Bioreactor example from bocop



function growth(s, mu2m, Ks)
    # MONOD
	return mu2m * s / (s+Ks)
end

function daynightgrowth(time, halfperiod, mubar)
    # light model: max^2 (0,sin) * mubar
	# DAY/NIGHT CYCLE: [0,2 halfperiod] rescaled to [0,2pi]
    pi = 3.141592653589793
	days = time / (halfperiod*2)
    tau = (days - floor(days)) * 2*pi
    mu = max(0,sin(tau))^2 * mubar
	return mu
end

#= METHANE PROBLEM
mu2 according to growth model, mu according to light model
time scale is [0,10] for 24h (day then night)=#
@def bioreactor_1 begin
    # constants
    beta = 1
    c = 2
    gamma = 1
    Ks = 0.05
    mu2m = 0.1
    mubar = 1
    r = 0.005
    T = 10

    t ∈ [ 0, T ], time
    x ∈ R³, state
    u ∈ R, control
    y = x[1]
    s = x[2]
    b = x[3]
    mu = daynightgrowth(t, T/2, mubar)
    mu2 = growth(s(t), mu2m, Ks)
    [0,0,0.001] ≤ x(t) ≤ [Inf, Inf, Inf]
    0 ≤ u(t) ≤ 1
    1 ≤ y(0) ≤ Inf
    1 ≤ b(0) ≤ Inf
    x(0)- x(T) == [0,0,0]
    ẋ(t) == [mu*y(t)/(1+y(t))-(r+u(t))*y(t),
            -mu2*b(t) + u(t)*beta*(gamma*y(t)-s(t)),
            (mu2-u(t)*beta)*b(t)]
    ∫(mu2*b(t)/(beta+c)) → max
end

#sol=solve(bioreactor_1)
#plot(sol)

return ((ocp=bioreactor_1, obj=0.614134, name="bioreactor_1day_periodic"))