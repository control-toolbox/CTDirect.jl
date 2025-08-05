using Plots

prob = truck_trailer()
sol = solve(prob.ocp; init=prob.init)

N = 50
x = zeros(N, 2)
for i in 1:N
    x[i, :] = state(sol)(i)[1:2]
end
plot(x[:, 2], x[:, 1])

# not obvious that the solution is the same as in the paper...
