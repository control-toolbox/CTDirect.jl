Jump / CTDirect comparison - algal bacterial problem
Ipopt version 3.14.17, running with linear solver MUMPS 5.7.3
Ipopt settings: tol=1e-8, mu_strategy=adaptive

Takeaways:
- Hessian seems to be handled differently by Jump, see nonzero values.
Maybe a less sparse but faster method is used ?
- convergence: iterations are different, total times are close.
- CTDirect still allocates way more memory (20 to 30x roughly, better than before).
For the biggest allocations, a significant time is passed during the AD phase, before Ipopt.
Likely some swap issue.

Note: update test CT with split time grid also, cf AD_backend branch


|                 | Trapeze (1000)           || Trapeze (5000)               ||
|                 | Jump   | CT      New     | Jump     | CT       | New      |
|-----------------|--------|--------|--------|----------|----------|----------|
|nnz jacobian     | 42006  | 42006  | 42006  | 210006   | 210006   | 210006   |
|nnz hessian      | 70000  | 12012  | 12012  | 350000   | 60012    | 60012    |
|variables        | 8008   | 8008   | 8008   | 40008    | 40008    | 40008    |
|lowerbound       | 6006   | 6006   | 6006   | 30006    | 30006    | 30006    |
|lower/upper      | 2002   | 2002   | 2002   | 10002    | 10002    | 10002    |
|equality         | 6006   | 6006   | 6006   | 30006    | 30006    | 30006    |
|-----------------|--------|--------|--------|----------|----------|----------|
|iterations       | 323    | 365    | 365    | 509      | 420      | 420      |
|objective        | 5.4522 | 5.4522 | 5.4522 | 5.4522   | 5.4522   | 5.4522   |
|structure        |        | noise  |        |          | noise    |          |
|allocations      | 336MB  | 4.5GB  | 4.7GB  | 3.0GB    | 49GB     | 50GB     |
|time (ipopt)     | 16(15) | 20(19) | 22(19) | 155(152) | 133(111) | 128(108) |


|                 | Gauss Legendre 2 (1000)  || Gauss Legendre 2 (5000)      ||
|                 | Jump   | CT     | New    | Jump     | CT       | New      |
|-----------------|--------|--------|--------|----------|----------|----------|
|nnz jacobian     |        | 138000 | 138000 |          | 690000   | 690000   |
|nnz hessian      |        | 81000  | 81000  |          | 405000   | 405000   |
|variables        |        | 20008  | 20008  |          | 100008   | 100008   |
|lowerbound       |        | 6006   | 6006   |          | 30006    | 30006    |
|lower/upper      |        | 2002   | 2002   |          | 10002    | 10002    |
|equality         |        | 18006  | 18006  |          | 90006    | 90006    |
|-----------------|--------|--------|--------|----------|----------|----------|
|iterations       |        | 100    | 100    |          | 111      | 111      |
|objective        |        | 5.4522 | 5.4522 |          | 5.4522   | 5.4522   |
|structure        |        | clean  |        |          | clean    |          |
|allocations      |        | 14.9GB | 14.5GB |          | 311GB    | 308GB    |
|time (ipopt)     |        | 41(29) | 39(31) |          | 393(143) | 371(140) |


