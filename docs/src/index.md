# CTDirect.jl

```@meta
CurrentModule =  CTDirect
```

The `CTDirect.jl` package is part of the [control-toolbox ecosystem](https://github.com/control-toolbox).

```@raw html
<img src="./assets/diagram.png" style="display: block; margin: 0 auto 20px auto;" width="320px">
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

Summary of the time discretization:

```math
\begin{array}{lcl}
t \in [t_0,t_f]   & \to & \{t_0, \ldots, t_N=t_f\}\\[0.2em]
x(\cdot),\, u(\cdot) & \to & X=\{x_0, \ldots, x_N, u_0, \ldots, u_N\} \\[1em]
\hline
\\
\text{criterion} & \to & \min\ g(x_0, x_N) \\[0.2em]
\text{dynamics}  & \to & x_{i+i} = x_i + ((t_{i+i} - t_i)/2)\, (f(t_i, x_i, u_i) + f(t_{i+1}, x_{i+1}, u_{i+1})) ~~ \text{(trapeze)}\\[0.2em]
\text{control constraints} &\to& \xi_l  \le  \xi(t_i, u_i)   \le \xi_u \\
\text{path constraints} &\to& \eta_l \le \eta(t_i, x_i)        \le \eta_u \\
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

We use packages from [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers) to solve the (NLP) problem.

!!! note " actual limitation"
    For the moment we have implemented
    - Only trapeze method for the discretization
    - Only Ipopt for the optimization soltware

!!! note "link with bocop"

    This package is equivalent to the [bocop](https://bocop.org) software 