Performance measurement
=======================

Simulation time
---------------

Due to historic reasons, SeisSol write in the stdout 3 times the
Simulation time: Wall time and CPU time are the same. Elapsed time (via
clock_gettime) takes time at a slightly earlier point, however the
difference should be minimal.

The "total time spent in compute kernels" is the total time accumulated
over all nodes without MPI communication. Basically this measures the
time in which all CPUs do some useful work.

FLOPS counter
-------------

SeisSol outputs the hardware (HW) and non-zero (NZ) flops. The first one
are the flops actually calculated by the CPU, and the second one is a
theoretical value which counts all flops that incorporate
multiplications or additions with zero. The reason why those are
distinct is that it is sometimes faster to compute a few extra zeros.
E.g. we do not always have efficient routines for sparse matrix matrix
multiplication, hence we use a routine for dense matrix matrix
multiplications which is faster but adds some extra work. So the
HW-GFLOP can be seen as machine utilization and the NZ-GFLOP can be seen
as actual work done.

Performance
-----------

One important value, which we usually publish in our papers, is the
performance in GFLOP/s/node. You can compute it with 
:math:`HW-(NZ-)GFLOP / #nodes / elapsed-time`.
You can compare this value with the publications in order to see if your
performance is ok.
