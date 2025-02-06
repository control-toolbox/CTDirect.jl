# Jump / CTDirect comparison - algal bacterial problem

Note that the problem is redefined for each method, jump and ctdirect.
Also, the Gauss Legendre 2 implementations use a piecewise constant control.

## Takeaways
- CTDirect still allocates at least x10 more memory, worsening for higher problem sizes
For the biggest allocations, a significant time is passed during the AD phase, before Ipopt. The runs marked with * spend half the time before optimization, likely some swap issue.
We note that Jump memory appears linear wrt steps for GL2, but a bit superlinear for Trapeze. CTDirect memory always increases superlinearly wrt steps.
- Hessian seems to be handled differently by Jump, see the higher nonzero values.
Maybe a less sparse but faster and less memory intensive method is used ? 
- convergence: objective and trajectory are similar, iterations differ, maybe due to the different hessian handling. Total computation times are similar for Trapeze and x2 to x5 slower for CTDirect for GL2, probably due to the memory effect.
- in terms of control structures, GL2 solutions are clean, Jump Trapeze solutions shows a bit of noise, while CTDirect Trapeze solutions are very noisy. How Jump manages to find a cleaner solution with Trapeze is unclear.

## Todo
- check on ipopt last iteration that tol is also 1e-8 for Jump
- test CTDirect with manual sparsity patterns
- investigate how jump finds a cleaner solution for trapeze discretization (print settings ?)

## Results: Jump vs CTDirect
See `test/jump_comparison.jl`
Ipopt details: `This is Ipopt version 3.14.14, running with linear solver MUMPS 5.6.2.`
Settings: tol=1e-8, mu_strategy=adaptive
```
Jump trapeze 1000:  15.889 s (7920529 allocations: 351.87 MiB)
Jump trapeze 2000:  59.236 s (23055275 allocations: 891.64 MiB)
Jump trapeze 5000:  128.857 s (55226845 allocations: 2.10 GiB)
Jump gauss_legendre_2 1000:  15.583 s (10998729 allocations: 726.48 MiB)
Jump gauss_legendre_2 2000:  26.836 s (21371405 allocations: 1.40 GiB)
Jump gauss_legendre_2 5000:  74.715 s (56343588 allocations: 3.58 GiB)
```

```
CTDirect (optimized) trapeze 1000:  19.976 s (46501061 allocations: 4.54 GiB)
CTDirect (optimized) trapeze 2000:  39.350 s (89302127 allocations: 12.26 GiB)
CTDirect (optimized) trapeze 5000:  127.653 s (267989402 allocations: 49.33 GiB)
CTDirect (optimized) gauss_legendre_2 1000:  30.309 s (36508333 allocations: 14.56 GiB)
CTDirect (optimized) gauss_legendre_2 2000:  90.069 s (85715676 allocations: 42.83 GiB)
CTDirect (optimized) gauss_legendre_2 5000:  293.751 s (159734254 allocations: 304.56 GiB)
```

```
CTDirect (manual) trapeze 1000:  49.673 s (66608733 allocations: 5.40 GiB)
CTDirect (manual) trapeze 2000:  112.605 s (140277493 allocations: 11.40 GiB)
CTDirect (manual) trapeze 5000:  336.893 s (418023859 allocations: 34.07 GiB)
CTDirect (manual) gauss_legendre_2 1000:  40.780 s (56179398 allocations: 4.20 GiB)
CTDirect (manual) gauss_legendre_2 2000:  94.828 s (127623066 allocations: 9.56 GiB)
CTDirect (manual) gauss_legendre_2 5000:  197.400 s (263619301 allocations: 19.84 GiB)
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
|allocations      | 352MB  | 4.5GB  | 5.50GB | 2.1GB    | 49GB     | 37GB     |
|time             | 17     | 20     | 50     | 126      | 136      | 354      |


## Details: Gauss Legendre 2 (1000 and 5000 steps)
+++redo ct

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
|allocations      | 726MB  | 14.6GB | 4.2GB  | 3.6GB    | 305GB    | 19.8GB   |
|time             | 15     | 28     | 40     | 77       | 291*     | 188      |

* half the time is before optimization, swap effect due to huge allocations ?
