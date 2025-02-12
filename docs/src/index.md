# CTDirect.jl

```@meta
CurrentModule =  CTDirect
```

The `CTDirect.jl` package is part of the [control-toolbox ecosystem](https://github.com/control-toolbox).

```mermaid
flowchart TD
O(<a href='https://control-toolbox.org/OptimalControl.jl/stable/'>OptimalControl</a>) --> B(<a href='https://control-toolbox.org/CTBase.jl/stable/'>CTBase</a>)
O --> D(<a href='https://control-toolbox.org/CTDirect.jl/stable/'>CTDirect</a>)
O --> F(<a href='https://control-toolbox.org/CTFlows.jl/stable/'>CTFlows</a>)
P(<a href='https://control-toolbox.org/CTProblems.jl/stable/'>CTProblems</a>) --> F
P --> B
F --> B
D --> B
style D fill:#FBF275
```

!!! note "Install"

    To install a package from the control-toolbox ecosystem, 
    please visit the [installation page](https://github.com/control-toolbox#installation).

An optimal control problem with fixed initial and final times, denoted (OCP), can be described as minimising the cost functional

```math
g(x(t_0), x(t_f)) + \int_{t_0}^{t_f} f^{0}(t, x(t), u(t))~\mathrm{d}t
```

where the state $x$ and the control $u$ are functions subject, for $t \in [t_0, t_f]$,
to the differential constraint

```math
   \dot{x}(t) = f(t, x(t), u(t))
```

and other constraints such as

```math
\begin{array}{llcll}
~\xi_l  &\le& \xi(t, u(t))        &\le& \xi_u, \\
\eta_l &\le& \eta(t, x(t))       &\le& \eta_u, \\
\psi_l &\le& \psi(t, x(t), u(t)) &\le& \psi_u, \\
\phi_l &\le& \phi(t_0, x(t_0), t_f, x(t_f)) &\le& \phi_u.
\end{array}
```

The so-called direct approach transforms the infinite dimensional optimal control problem (OCP) into a finite dimensional optimization problem (NLP). This is done by a discretization in time by Runge-Kutta methods applied to the state and control variables, as well as the dynamics equation. These methods are usually less precise than indirect methods based on [Pontryagin’s Maximum Principle](https://en.wikipedia.org/w/index.php?title=Pontryagin's_maximum_principle&oldid=1160355192), but more robust with respect to the initialization. Also, they are more straightforward to apply, hence their wide use in industrial applications. We refer the reader to for instance[^1] and [^2] for more details on direct transcription methods and NLP algorithms.

[^1]: J. T. Betts. Practical methods for optimal control using nonlinear programming. Society for Industrial and Applied Mathematics (SIAM), Philadelphia, PA, 2001.

[^2]: J. Nocedal and S.J. Wright. Numerical optimization. Springer-Verlag, New York, 1999.****

Example of the time discretization by the trapezoidal rule:

```math
\begin{array}{lcl}
t \in [t_0,t_f]   & \to & \{t_0, \ldots, t_N=t_f\}\\[0.2em]
x(\cdot),\, u(\cdot) & \to & X=\{x_0, \ldots, x_N, u_0, \ldots, u_N\} \\[1em]
\hline
\\
\text{step} & \to & h = (t_f-t_0)/N\\[0.2em]
\text{criterion} & \to & \min\ g(x_0, x_N) \\[0.2em]
\text{dynamics}  & \to & x_{i+i} = x_i + (h/2)\, (f(t_i, x_i, u_i) + f(t_{i+1}, x_{i+1}, u_{i+1})) \\[0.2em]
\text{control constraints} &\to& \xi_l  \le  \xi(t_i, u_i)   \le \xi_u \\[0.2em]
\text{path constraints} &\to& \eta_l \le \eta(t_i, x_i)        \le \eta_u \\[0.2em]
\text{mixed constraints} &\to& \psi_l \le \psi(t_i, x_i, u_i) \le \psi_u \\[0.2em]
\text{limit conditions} &\to& \phi_l \le \phi(x_0, x_N) \le \phi_u
\end{array}
```

We therefore obtain a nonlinear programming problem on the discretized state and control variables of the general form:

```math
(NLP)\quad \left\{
\begin{array}{lr}
\min \ F(X) \\
LB \le C(X) \le UB
\end{array}
\right.
```

Solving the (NLP) problem is done using packages from [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers), with Ipopt as the default solver.

On the input side of this package, we use an [`OptimalControlModel`](@ref) structure from CTBase to define the (OCP).

The direct transcription to build the (NLP) can use discretization schemes such as trapeze (default), midpoint, or Gauss-Legendre collocations.

!!! note "Related packages"

    This package is equivalent to the [bocop](https://www.bocop.org) software.
