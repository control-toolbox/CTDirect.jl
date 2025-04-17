# Jump / CTDirect comparison - algal bacterial problem

Note that the problem is redefined for each method, jump and ctdirect.
Also, the Gauss Legendre 2 implementations use a piecewise constant control.

## Takeaways
- convergence: objective and trajectory are similar, iterations differ, maybe due to the different hessian handling. Total computation times for Trapeze are similar for Jump and CTDirect (optimized backend). For GL2, Jump is faster than both CTDirect versions.
- CTDirect with manual sparsity patterns still allocates 5 to 10 times more memory than Jump, but scales better than the automatic sparsity detection. Manual mode becomes faster than automatic mode for GL2 above 2000 steps, while being slower for Trapeze even at 5000 steps.
- Hessian is handled differently by Jump, with a less sparse but faster method.
- in terms of control structures, GL2 solutions are clean, Jump Trapeze solutions shows a bit of noise, while CTDirect Trapeze solutions are very noisy. How Jump manages to find a cleaner solution with Trapeze is unclear.

## Todo
- investigate how jump finds a cleaner solution for trapeze discretization (print settings ?)

## Results: Jump vs CTDirect
See `test/jump_comparison.jl`
Ipopt details: `This is Ipopt version 3.14.14, running with linear solver MUMPS 5.6.2.`
Settings: tol=1e-8, mu_strategy=adaptive
```
Jump trapeze 1000:  15.633 s (7920529 allocations: 351.87 MiB)
Jump trapeze 2000:  57.472 s (23055275 allocations: 891.64 MiB)
Jump trapeze 5000:  124.108 s (55226845 allocations: 2.10 GiB)
Jump gauss_legendre_2 1000:  15.250 s (10998729 allocations: 726.48 MiB)
Jump gauss_legendre_2 2000:  26.686 s (21371405 allocations: 1.40 GiB)
Jump gauss_legendre_2 5000:  76.455 s (56343588 allocations: 3.58 GiB)

Jump trapeze 1000:  15.960 s (9142234 allocations: 636.63 MiB)
Jump trapeze 2000:  57.480 s (27887980 allocations: 1.86 GiB)
Jump trapeze 5000:  125.863 s (66492550 allocations: 4.42 GiB)
Jump gauss_legendre_2 1000:  15.883 s (11302329 allocations: 936.87 MiB)
Jump gauss_legendre_2 2000:  27.211 s (21755005 allocations: 1.77 GiB)
Jump gauss_legendre_2 5000:  78.946 s (58251188 allocations: 4.66 GiB)


CTDirect (optimized) trapeze 1000:  17.941 s (45041839 allocations: 4.47 GiB)
CTDirect (optimized) trapeze 2000:  35.811 s (86384903 allocations: 12.12 GiB)
CTDirect (optimized) trapeze 5000:  124.451 s (260698176 allocations: 48.99 GiB)
CTDirect (optimized) gauss_legendre_2 1000:  25.087 s (32172053 allocations: 14.36 GiB)
CTDirect (optimized) gauss_legendre_2 2000:  76.272 s (77043394 allocations: 42.43 GiB)
CTDirect (optimized) gauss_legendre_2 5000:  281.000 s (138053972 allocations: 303.56 GiB)

CTDirect (optimized) trapeze 1000:  14.642 s (5846082 allocations: 2.71 GiB)
CTDirect (optimized) trapeze 2000:  36.932 s (13040418 allocations: 9.19 GiB)
CTDirect (optimized) trapeze 5000:  121.042 s (38176659 allocations: 40.05 GiB)
CTDirect (optimized) gauss_legendre_2 1000:  30.150 s (9681285 allocations: 13.66 GiB)
CTDirect (optimized) gauss_legendre_2 2000:  67.987 s (17028256 allocations: 39.65 GiB)
CTDirect (optimized) gauss_legendre_2 5000:

CTDirect (manual) trapeze 1000:  47.442 s (65149511 allocations: 5.33 GiB)
CTDirect (manual) trapeze 2000:  105.765 s (137360269 allocations: 11.26 GiB)
CTDirect (manual) trapeze 5000:  324.819 s (410732633 allocations: 33.73 GiB)
CTDirect (manual) gauss_legendre_2 1000:  38.147 s (51843118 allocations: 4.00 GiB)
CTDirect (manual) gauss_legendre_2 2000:  87.863 s (118950784 allocations: 9.16 GiB)
CTDirect (manual) gauss_legendre_2 5000:  187.292 s (241939016 allocations: 18.85 GiB)

CTDirect (manual) trapeze 1000:  48.345 s (8575383 allocations: 3.10 GiB)
CTDirect (manual) trapeze 2000:  106.494 s (17172721 allocations: 6.23 GiB)
CTDirect (manual) trapeze 5000:  317.398 s (48925210 allocations: 17.97 GiB)
CTDirect (manual) gauss_legendre_2 1000:  49.515 s (10792741 allocations: 2.73 GiB)
CTDirect (manual) gauss_legendre_2 2000:  58.446 s (14574145 allocations: 3.37 GiB)
CTDirect (manual) gauss_legendre_2 5000:  236.905 s (50702300 allocations: 12.47 GiB)

```

## Details: Trapeze (1000 and 5000 steps)

|                 | Jump   | CT     | Manual | Jump     | CT       | Manual   |
|-----------------|--------|--------|--------|----------|----------|----------|
|nnz jacobian     | 42006  | 42006  | 96076  | 210006   | 210006   | 480072   |
|nnz hessian      | 74000  | 12012  | 100072 | 370000   | 60012    | 500072   |
|variables        | 8008   | 8008   | 8008   | 40008    | 40008    | 40008    |
|lowerbound       | 6006   | 6006   | 6006   | 30006    | 30006    | 30006    |
|lower/upper      | 2002   | 2002   | 2002   | 10002    | 10002    | 10002    |
|equality         | 6006   | 6006   | 6006   | 30006    | 30006    | 30006    |
|iterations       | 334    | 365    | 333    | 517      | 420      | 419      |
|objective        | 5.4522 | 5.4522 | 5.4522 | 5.4522   | 5.4522   | 5.4522   |
|structure        | ok     | noisy  | noisy  | ok       | noisy    | noisy    |

## Details: Gauss Legendre 2 (1000 and 5000 steps)

|                 | Jump   | CT     | Manual | Jump     | CT       | Manual   |
|-----------------|--------|--------|--------|----------|----------|----------|
|nnz jacobian     | 118006 | 118006 | 384072 | 590006   | 590006   | 1920072  |
|nnz hessian      | 322000 | 63000  | 210072 | 1610000  | 315000   | 1050072  |
|variables        | 20006  | 20008  | 20008  | 100006   | 100008   | 100008   |
|lowerbound       | 6006   | 6006   | 6006   | 3006     | 30006    | 30006    |
|lower/upper      | 2000   | 2002   | 2002   | 10000    | 10002    | 10002    |
|equality         | 18006  | 18006  | 18006  | 90006    | 90006    | 90006    |
|iterations       | 117    | 95     | 93     | 146      | 78       | 86       |
|objective        | 5.4522 | 5.4522 | 5.4522 | 5.4522   | 5.4522   | 5.4522   |
|structure        | clean  | clean  | clean  | clean    | clean    | clean    |

* half the time is before optimization, swap effect due to huge allocations ?
