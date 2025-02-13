# NLP solvers comparison

Ipopt: mumps linear solver
MadNLP: umfpack linear solver
Knitro: 14.2.0, default options except feas/opt tolerances set to CTDirect default value

Benchmarks



Trapeze

```
julia> bench(;nlp_solver=:ipopt)
Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
Grid size list: [250, 500, 1000, 2500, 5000]
Grid size    250: time (s) =    0.6
Grid size    500: time (s) =    1.5
Grid size   1000: time (s) =    4.5
Grid size   2500: time (s) =   18.4
Grid size   5000: time (s) =   69.7

julia> bench(;nlp_solver=:madnlp)
Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
Grid size list: [250, 500, 1000, 2500, 5000]
Grid size    250: time (s) =    0.8
Grid size    500: time (s) =    1.9
Grid size   1000: time (s) =    5.3
Grid size   2500: time (s) =   22.8
Grid size   5000: time (s) =   82.0

julia> bench(;nlp_solver=:knitro)
Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
Grid size list: [250, 500, 1000, 2500, 5000]
Grid size    250: time (s) =    0.6
Grid size    500: time (s) =    1.3
Grid size   1000: time (s) =    3.6
Grid size   2500: time (s) =   16.1
Grid size   5000: time (s) =   58.5
```

Gauss Legendre 2 with manual sparsity patterns
```
julia> bench(;disc_method=:gauss_legendre_2, adnlp_backend=:manual)
Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
Grid size list: [250, 500, 1000, 2500, 5000]
Grid size    250: time (s) =    3.5
Grid size    500: time (s) =    9.5
Grid size   1000: time (s) =   18.7
Grid size   2500: time (s) =   57.2
Grid size   5000: time (s) =  122.2

julia> bench(;disc_method=:gauss_legendre_2, adnlp_backend=:manual, nlp_solver=:madnlp)
Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
Grid size list: [250, 500, 1000, 2500, 5000]
Grid size    250: time (s) =   15.9
Grid size    500: time (s) =   37.5
Grid size   1000: time (s) =   83.6
Grid size   2500: time (s) =  759.3

julia> bench(;disc_method=:gauss_legendre_2, adnlp_backend=:manual, nlp_solver=:knitro)
Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
Grid size list: [250, 500, 1000, 2500, 5000]
Grid size    250: time (s) =   12.0
Grid size    500: time (s) =   10.8
Grid size   1000: time (s) =   47.4
Grid size   2500: time (s) =  132.1
Grid size   5000: time (s) =  341.7
```