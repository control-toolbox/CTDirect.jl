# Benchmark for different AD backends
The backend for ADNLPModels can be set in transcription / solve calls with the option `adnlp_backend = ...`. Possible values are, including the predefined backups for ADNLPModels (*) :
- :optimized* (default for CTDirect) : Forward for Jacobian, Reverse for Gradient and Hessian
- :default* : Forward for everything (much slower)
- :manual : give to ADNLPModels the sparse pattern for Jacobian and Hessian (use same Forward / Reverse settings as the :optimized predefined backend)  
- :enzyme* : Enzyme (not working)
- :zygote* : Zygote (not working)

## Errors:
- enzyme gives correct nonzero counts for Jacobian and Hessian, but fails with
```ERROR: Constant memory is stored (or returned) to a differentiable variable.
As a result, Enzyme cannot provably ensure correctness and throws this error.
This might be due to the use of a constant variable as temporary storage for active memory (https://enzyme.mit.edu/julia/stable/faq/#Runtime-Activity).
If Enzyme should be able to prove this use non-differentable, open an issue!
To work around this issue, either:
 a) rewrite this variable to not be conditionally active (fastest, but requires a code change), or
 b) set the Enzyme mode to turn on runtime activity (e.g. autodiff(set_runtime_activity(Reverse), ...) ). This will maintain correctness, but may slightly reduce performance.```
 Error apparently occurs when calling the boundary conditions.```
- zygote gives incorrect (huge) nonzero counts then also fails with an error message. 

## Tests:
Using CTDirect benchmark function bench()
Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "robbins", "simple_integrator", "vanderpol"]

Takeaways:
- the optimized backend (with ReverseDiff for Hessian) is much better than full ForwardDiff.
- manual sparse pattern seems to give better performance for larger problems, likely because of the hugely increasing cost of computing the Hessian sparsity in terms of allocations and time (cf also comparison with Jump that uses a different, less sparse but faster Hessian).

| Trapeze | default | optimized | manual  |
|---------|---------|-----------|---------|
| 250     | 49.7    | 0.9       | 1.5     |
| 500     |         | 2.4       | 3.4     |
| 1000    |         | 6.2       | 6.4     |
| 2500    |         | 24.7      | 20.6    |
| 5000    |         |           | 50.0    |
| 7500    |         |           | 61.2    |
| 10000   |         |           |    |


Sparsity details: goddard_all Trapeze (1000 and 10000 steps)

| transcription | optimized | manual  | optimized | manual |
|---------------|-----------|---------|-----------|--------|
| NLP vars      | 4005      | 4005    | 40005     | 40005  |
| NLP cons      | 6007      | 6007    | 60007     | 60007  |
| Hess nnz      | 11011     | 30024   | 110011    | 300024 |
| H sparsity    | 99.86%    | 99.63%  | 99.99%    | 99.96% |
| Jac nnz       | 28011     | 42043   | 280011    | 420043 |
| J sparsity    | 99.88%    | 99.83%  | 99.99%    | 99.98% |
| allocs        | 1.16GB    | 106MB   | 71.56GB   | 4.55GB |
| time          | 750ms     | 85ms    | 64.7s**   | 3.8s   |
|---------------|-----------|---------|-----------|--------|
| solve         | optimized | manual  | optimized | manual |
| iterations    | 42        | 28      | 51        | 29     |
| allocs        | 2.0GB     | 1.2GB   | 87.5GB    | 16.9GB |
| time          | 2.5s      | 2.5s    | 151.0s*** | 42.4s  |

** hessian accounts for 59 out of total 65s
*** building the hessian is one third of the total solve time !

## Todo:
- try to disable some unused (?) parts such as hprod ? (according to show_time info the impact may be small)
- add pattern structure for midpoint and IRK schemes
- redo tests on algal_bacterial problem, including Jump
- add some tests for different backends in test_misc
- reuse ADNLPModels functions to get block sparsity patterns then rebuild full patterns ?

