Benchmark CTDirect for different AD backends (keyword adnlp_backend = ...)
- :default : ForwardDiff
- :optimized (default for CTDirect) : ForwardDiff for Jacobian / ReverseDiff for Hessian
- :enzyme : Enzyme
- :zygote : Zygote

Using CTDirect benchmark function bench()
Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "robbins", "simple_integrator", "vanderpol"]

Takeaways:
- optimized backend (with ReverseDiff for Hessian) is much better than full ForwardDiff.
- enzyme gives correct nonzero counts for Jacobian and Hessian, but fails with
```ERROR: Constant memory is stored (or returned) to a differentiable variable.
As a result, Enzyme cannot provably ensure correctness and throws this error.
This might be due to the use of a constant variable as temporary storage for active memory (https://enzyme.mit.edu/julia/stable/faq/#Runtime-Activity).
If Enzyme should be able to prove this use non-differentable, open an issue!
To work around this issue, either:
 a) rewrite this variable to not be conditionally active (fastest, but requires a code change), or
 b) set the Enzyme mode to turn on runtime activity (e.g. autodiff(set_runtime_activity(Reverse), ...) ). This will maintain correctness, but may slightly reduce performance.```
 Error apparently occurs when calling the boundary conditions.
- zygote gives incorrect (huge) nonzero counts then also fails with an error message. 


| Trapeze | default | optimized |
|---------|---------|-----------|
| 250     | 43.3    | 1.0       |
| 500     | 176.2   | 2.4       |
| 1000    | 926.0   | 7.1       |
| 2500    |         | 31.8      |
| 5000    |         |           |

Midpoint

GaussLegendre2
