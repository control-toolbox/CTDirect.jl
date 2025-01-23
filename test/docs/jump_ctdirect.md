Jump / CTDirect comparison - algal bacterial problem (see test/jump_comparison.jl)
Ipopt version 3.14.17, running with linear solver MUMPS 5.7.3
Ipopt settings: tol=1e-8, mu_strategy=adaptive
Note that the problem is redefined for each method: jump, ctdirect and ctdirect new model.

Also, the Gauss Legendre 2 implementations for Jump and CTDirect handle the control slightly differently: jump version uses a piecewise constant control while CTDirect uses a piecewise linear control. So the discretization is not exactly identical, and the number of nonzero elements in derivatives will change, see **. This also means that the CTDirect version has mN+1 controls instead of N for the Jump version. Neither versions use the formulation where stage controls are directly treated as variables, however for a 2-stage method the number of control variables will be similar.

Takeaways:
- CTDirect still allocates way more memory (roughly x15, better than before).
For the biggest allocations, a significant time is passed during the AD phase, before Ipopt. The runs marked with * spend half the time before optimization, likely some swap issue.
- Hessian seems to be handled differently by Jump, see the higher nonzero values.
Maybe a less sparse but faster and less memory intensive method is used ? 
- convergence: iterations are different, maybe due to the different hessian handling.
Total computation times are close.

Todo:
- disable Hessian (use ipopt finite differences) and compare memory allocations and convergence
- piecewise constant control for GL2 in CTDirect




|                 | Trapeze (1000)           || Trapeze (5000)               ||
|                 | Jump   | CT      New     | Jump     | CT       | New      |
|-----------------|--------|--------|--------|----------|----------|----------|
|nnz jacobian     | 42006  | 42006  | 42006  | 210006   | 210006   | 210006   |
|nnz hessian      | 70000  | 12012  | 12012  | 350000   | 60012    | 60012    |
|variables        | 8008   | 8008   | 8008   | 40008    | 40008    | 40008    |
|lowerbound       | 6006   | 6006   | 6006   | 30006    | 30006    | 30006    |
|lower/upper      | 2002   | 2002   | 2002   | 10002    | 10002    | 10002    |
|equality         | 6006   | 6006   | 6006   | 30006    | 30006    | 30006    |
|iterations       | 323    | 365    | 365    | 509      | 420      | 420      |
|objective        | 5.4522 | 5.4522 | 5.4522 | 5.4522   | 5.4522   | 5.4522   |
|structure        |        | noisy  |        |          | noisy    |          |
|-----------------|--------|--------|--------|----------|----------|----------|
|allocations      | 340MB  | 4.5GB  | 4.7GB  | 3.0GB    | 49GB     | 50GB     |
|time             | 19     | 20     | 22     | 155      | 133      | 128      |


|                 | Gauss Legendre 2 (1000)  || Gauss Legendre 2 (5000)      ||
|                 | Jump   | CT     | New    | Jump     | CT       | New      |
|-----------------|--------|--------|--------|----------|----------|----------|
|nnz jacobian     | 118006 | 138000** | 138000** |          | 690000**   | 690000**   |
|nnz hessian      | 322000 | 81000  | 81000  |          | 405000   | 405000   |
|variables        | 20006  | 20008  | 20008  |          | 100008   | 100008   |
|lowerbound       | 6006   | 6006   | 6006   |          | 30006    | 30006    |
|lower/upper      | 2000   | 2002   | 2002   |          | 10002    | 10002    |
|equality         | 18006  | 18006  | 18006  |          | 90006    | 90006    |
|iterations       |        | 100    | 100    |          | 111      | 111      |
|objective        |        | 5.4522 | 5.4522 |          | 5.4522   | 5.4522   |
|structure        |        | clean  |        |          | clean    |          |
|-----------------|--------|--------|--------|----------|----------|----------|
|allocations      |        | 14.9GB | 14.5GB |          | 313GB    | 308GB    |
|time             |        | 35     | 39     |          | 382*     | 371*     |


* swap effect due to huge allocations ?
** slightly different implementation for the stage controls


Jump trapeze 1000:  17.840 s (6779482 allocations: 340.27 MiB)
Jump trapeze 2000:  58.605 s (18798587 allocations: 917.26 MiB)
Jump trapeze 5000:  149.540 s (47597571 allocations: 3.00 GiB)
Jump gauss_legendre_2 1000:  15.528 s (10988856 allocations: 726.32 MiB)
Jump gauss_legendre_2 2000:  26.995 s (21345532 allocations: 1.40 GiB)
Jump gauss_legendre_2 5000:  75.047 s (56269715 allocations: 3.57 GiB)
CTDirect trapeze 1000:  20.110 s (46501059 allocations: 4.54 GiB)
CTDirect trapeze 2000:  41.097 s (89302125 allocations: 12.26 GiB)
CTDirect trapeze 5000:  133.268 s (267989400 allocations: 49.33 GiB)
piecewise constant control

piecewise linear control
CTDirect gauss_legendre_2 1000:  37.220 s (39259673 allocations: 15.12 GiB)
CTDirect gauss_legendre_2 2000:  112.687 s (104950745 allocations: 45.37 GiB)
CTDirect gauss_legendre_2 5000:  363.224 s (211848763 allocations: 313.02 GiB)