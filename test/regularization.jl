# test regularization (tychonov)
using Plots

# algal_bacterial
# nb continuation does not seem really useful (slower)
# Q. how to adjust epsilon automatically ? use some relative weight ?
# ie int eps u² instead of eps int u² ? meh...
sol0 = solve(algal_bacterial().ocp)
plot(sol0)
sol = solve(algal_bacterial(1e-2).ocp)
plot!(sol)
