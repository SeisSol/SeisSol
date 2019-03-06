Computing time vs order of accuracy
===================================

Using a higher order implies using smaller time steps and using more
basis functions.
The expected increase in computed time can then be estimated as follow:

| Decrease of the time steps width by a factor: :math:`(2N+1)/(2n+1)`
| Increase of the number of basis functions by a factor:
  :math:`(N+1)(N+2)(N+3)/((n+1)(n+2)(n+3))`
| with:
| n: initial order of accuracy -1
| N: ne w order of accuracy -1

Here is an example of the theoretical increase of the compute time
relative to order 3

========================================= = === === ===== ===== == ==== ====
\                                         3 4   5   6     7     8  9    10
========================================= = === === ===== ===== == ==== ====
Increase due to time steps                1 1.4 1.8 2.2   2.6   3  3.4  3.8
Increase due to number of basis functions 1 2   3.5 5.6   8.4   12 16.5 22
theoretical increase relative to order 3  1 2.8 6.3 12.32 21.84 36 56.1 83.6
observed increase (on SM2)                1 2.0 4.6 11.5               
========================================= = === === ===== ===== == ==== ====

The last line show the observed time increase on SuperMUC Phase 2 on a small
run with dynamic rupture and LTS-DR. Low order calculation (up to order
5 included) are memory bounds, and are then less efficient. As a
consequence, the higher-order simulations cost less than expected in
comparison with order 3.
