# Benchmark CTDirect for different AD backends
The backend for ADNLPModels can be set in transcription / solve calls with the option `adnlp_backend = ...`
- :optimized (default for CTDirect) : ForwardDiff for Jacobian / ReverseDiff for Hessian
- :default : ForwardDiff (much slower)
- :manual : sparse pattern for Jacobian / Hessian given to ADNLPModels 
- :enzyme : Enzyme
- :zygote : Zygote

## Errors:
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

## Tests:
Using CTDirect benchmark function bench()
Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "robbins", "simple_integrator", "vanderpol"]

Takeaways:
- optimized backend (with ReverseDiff for Hessian) is much better than full ForwardDiff.
- manual sparse pattern seems to give better performance for larger problems, likely because of the hugely increasing cost of computing the Hessian sparsity in terms of allocations and time (cf also comparison with Jump that uses a different, less sparse but faster Hessian).

| Trapeze | default | optimized | manual  |
|---------|---------|-----------|---------|
| 250     | 43.3    | 1.0       | 1.6     |
| 500     | 176.2   | 2.6       | 3.8     |
| 1000    | 926.0   | 7.3       | 7.7     |
| 2500    |         | 29.2      | 30.5    |
| 5000    |         | 108.5     | 84.1    |
| 6000    |         |           | 94.2    |
| 7000    |         |           | 130.8   |
| 8000    |         |           | 166.0   |
| 9000    |         |           | 201.7   |
| 10000   |         | 1252.4    | 154.0   |


Sparsity details: goddard_all Trapeze (1000)
| transcription | :optimized | :manual |
|---------------|------------|---------|
| NLP vars      | 4005       | 4005    |
| NLP cons      | 6007       | 6007    |
| Hess nnz      | 11011      | 30024   |
| H sparsity    | 99.86%     | 99.63%  |
| Jac nnz       | 28011      | 42043   |
| J sparsity    | 99.88%     | 99.83%  |
| allocs        | 1.16GB     | 106MB   |
| time          | 750ms      | 85ms    |
|---------------|------------|---------|
| solve         | :optimized | :manual |
| iterations    | 42         | 28      |
| allocs        | 2.0GB      | 1.2GB   |
| time          | 2.5s       | 2.5s    |

Sparsity details: algal_bacterial Trapeze (1000)
| transcription | :optimized | :manual | Jump
|---------------|------------|---------|
| NLP vars      |            |         |
| NLP cons      |            |         |
| Hess nnz      |            |         |
| H sparsity    |            |         |
| Jac nnz       |            |         |
| J sparsity    |            |         |
| allocs        |            |         |
| time          |            |         |
|---------------|------------|---------|
| solve         | :optimized | :manual |
| iterations    |            |         |
| allocs        |            |         |
| time          |            |         |

## Todo:
- check all specific backends (jprod etc) and set them as in :optimized
- add pattern structure for midpoint and IRK schemes
- redo tests on algal_bacterial problem, including Jump
- reuse ADNLPModels functions to get block sparsity patterns then rebuild full patterns ?

