# NLP solvers comparison

Solvers tested:
- Ipopt: 3.14.14, linear solver MUMPS 5.6.2 (option mu_strategy="adaptive")
- MadNLP: v0.8.5, linear solver umfpack
- Knitro: 14.2.0, using the Interior-Point/Barrier Direct algorithm.
All solvers use the same primal / dual tolerances, set to CTDirect default.

Benchmarks for trapeze and Gauss Legendre 2 discretizations
- Problem list: ["beam", "double_integrator_mintf", "double_integrator_minenergy", "double_integrator_freet0tf", "fuller", "goddard", "goddard_all", "jackson", "simple_integrator", "vanderpol"]
- Grid size list: [250, 500, 1000, 2500, 5000]

Trapeze: performance is close for the 3 solvers
```
julia> bench(;nlp_solver=:ipopt)
Grid size    250: time (s) =    0.6
Grid size    500: time (s) =    1.4
Grid size   1000: time (s) =    3.8
Grid size   2500: time (s) =   18.0
Grid size   5000: time (s) =   74.8

julia> bench(;nlp_solver=:madnlp)
Grid size    250: time (s) =    0.7
Grid size    500: time (s) =    1.7
Grid size   1000: time (s) =    4.7
Grid size   2500: time (s) =   19.4
Grid size   5000: time (s) =   74.1

julia> bench(;nlp_solver=:knitro)
Grid size    250: time (s) =    0.5
Grid size    500: time (s) =    1.3
Grid size   1000: time (s) =    3.6
Grid size   2500: time (s) =   16.0
Grid size   5000: time (s) =   58.2
```

Gauss Legendre 2 with manual sparsity patterns
```
julia> bench(;disc_method=:gauss_legendre_2, adnlp_backend=:manual, nlp_solver=:ipopt)
Grid size    250: time (s) =    3.2
Grid size    500: time (s) =    7.8
Grid size   1000: time (s) =   18.5
Grid size   2500: time (s) =   56.4
Grid size   5000: time (s) =  130.8

julia> bench(;disc_method=:gauss_legendre_2, adnlp_backend=:manual, nlp_solver=:madnlp)
Grid size    250: time (s) =   14.3
Grid size    500: time (s) =   36.6
Grid size   1000: time (s) =   85.9
Grid size   2500: time (s) =  832.7

julia> bench(;disc_method=:gauss_legendre_2, adnlp_backend=:manual, nlp_solver=:knitro)
Grid size    250: time (s) =    6.8
Grid size    500: time (s) =   10.5
Grid size   1000: time (s) =   44.4
Grid size   2500: time (s) =  120.3
Grid size   5000: time (s) =  338.4
```