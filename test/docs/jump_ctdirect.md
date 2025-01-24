# Jump / CTDirect comparison - algal bacterial problem

Note that the problem is redefined for each method: jump, ctdirect and ctdirect new model.
Also, the Gauss Legendre 2 implementations for Jump and CTDirect here use a piecewise constant control (default for CTDirect would have been piecewise linear).

## Takeaways
- CTDirect still allocates at least x10 more memory, worsening for higher problem sizes
For the biggest allocations, a significant time is passed during the AD phase, before Ipopt. The runs marked with * spend half the time before optimization, likely some swap issue.
We note that Jump memory appears linear wrt steps for GL2, but a bit superlinear for Trapeze. CTDirect memory always increases superlinearly wrt steps.
- Hessian seems to be handled differently by Jump, see the higher nonzero values.
Maybe a less sparse but faster and less memory intensive method is used ? 
- convergence: objective and trajectory are similar, iterations differ, maybe due to the different hessian handling. Total computation times are similar for Trapeze and x2 to x5 slower for CTDirect for GL2, probably due to the memory effect.
- for GL2, Jump and CTDirect have slightly different nonzero counts for the Jacobian
- in terms of control structures, GL2 solutions are clean, Jump Trapeze solutions shows a bit of noise, while CTDirect Trapeze solutions are very noisy.

## Todo
- can we have linear memory wrt steps for Jump / Trapeze ?
- disable Hessian (in AD model then use ipopt limited memory option ?) and compare memory allocations and convergence. Find more details on the Hessian in Jump.

## Results: Jump vs CTDirect
See `test/jump_comparison.jl`
Ipopt details: `Ipopt version 3.14.17, running with linear solver MUMPS 5.7.3`
Settings: tol=1e-8, mu_strategy=adaptive

```
Jump trapeze 1000:  17.029 s (7920527 allocations: 351.87 MiB)
Jump trapeze 2000:  56.928 s (23055273 allocations: 891.64 MiB)
Jump trapeze 5000:  125.798 s (55226843 allocations: 2.10 GiB)
Jump gauss_legendre_2 1000:  15.398 s (10988856 allocations: 726.32 MiB)
Jump gauss_legendre_2 2000:  27.401 s (21345532 allocations: 1.40 GiB)
Jump gauss_legendre_2 5000:  76.593 s (56269715 allocations: 3.57 GiB)
```
Trapeze
```
CTDirect trapeze 1000:  20.110 s (46501059 allocations: 4.54 GiB)
CTDirect trapeze 2000:  41.097 s (89302125 allocations: 12.26 GiB)
CTDirect trapeze 5000:  133.268 s (267989400 allocations: 49.33 GiB)
```
new model
```
CTDirect trapeze 1000:  18.833 s (3227130 allocations: 2.27 GiB)
CTDirect trapeze 2000:  38.233 s (6469103 allocations: 7.90 GiB)
CTDirect trapeze 5000:  119.554 s (16850293 allocations: 36.15 GiB)
```

GL2 piecewise constant control
```
CTDirect gauss_legendre_2 1000:  33.181 s (37843213 allocations: 14.79 GiB)
CTDirect gauss_legendre_2 2000:  82.605 s (82766476 allocations: 43.19 GiB)
CTDirect gauss_legendre_2 5000:  356.338 s (221161426 allocations: 312.28 GiB)
```
new model
```
CTDirect gauss_legendre_2 1000:  34.719 s (13521528 allocations: 13.47 GiB)
CTDirect gauss_legendre_2 2000:  81.367 s (28519019 allocations: 40.26 GiB)
CTDirect gauss_legendre_2 5000:  367.255 s (72472563 allocations: 304.25 GiB)
```

GL2 piecewise linear control
```
CTDirect gauss_legendre_2 1000:  37.220 s (39259673 allocations: 15.12 GiB)
CTDirect gauss_legendre_2 2000:  112.687 s (104950745 allocations: 45.37 GiB)
CTDirect gauss_legendre_2 5000:  363.224 s (211848763 allocations: 313.02 GiB)
```

## Details: Trapeze (1000 and 5000 steps)

|                 | Jump   | CT     | New    | Jump     | CT       | New      |
|-----------------|--------|--------|--------|----------|----------|----------|
|nnz jacobian     | 42006  | 42006  | 42006  | 210006   | 210006   | 210006   |
|nnz hessian      | 74000  | 12012  | 12012  | 370000   | 60012    | 60012    |
|variables        | 8008   | 8008   | 8008   | 40008    | 40008    | 40008    |
|lowerbound       | 6006   | 6006   | 6006   | 30006    | 30006    | 30006    |
|lower/upper      | 2002   | 2002   | 2002   | 10002    | 10002    | 10002    |
|equality         | 6006   | 6006   | 6006   | 30006    | 30006    | 30006    |
|iterations       | 334    | 365    | 365    | 517      | 420      | 420      |
|objective        | 5.4522 | 5.4522 | 5.4522 | 5.4522   | 5.4522   | 5.4522   |
|structure        | ok     | noisy  |        | ok       | noisy    |          |
|allocations      | 352MB  | 4.5GB  |        | 2.1GB    | 49GB     |          |
|time             | 17     | 20     |        | 126      | 136      |          |


## Details: Gauss Legendre 2 (1000 and 5000 steps)

|                 | Jump   | CT     | New    | Jump     | CT       | New      |
|-----------------|--------|--------|--------|----------|----------|----------|
|nnz jacobian     | 118006 | 124000 |        | 590006   | 620000   |          |
|nnz hessian      | 322000 | 63000  |        | 1610000  | 315000   |          |
|variables        | 20006  | 20008  |   | 100006   | 100008   |    |
|lowerbound       | 6006   | 6006   |    | 3006     | 30006    |     |
|lower/upper      | 2000   | 2002   |    | 10000    | 10002    |     |
|equality         | 18006  | 18006  |   | 90006    | 90006    |     |
|iterations       | 117    | 96     |        | 146      | 119      |          |
|objective        | 5.4522 | 5.4522 |  | 5.4522   | 5.4522   |    |
|structure        | clean  | clean  |        | clean    | clean    |          |
|allocations      | 726MB  | 14.8GB |        | 3.6GB    | 312GB    |          |
|time             | 15     | 33     |        | 77       | 356*     |          |

* half the time is before optimization, swap effect due to huge allocations ?



