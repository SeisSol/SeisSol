..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Performance measurement
=======================

Simulation time
---------------

SeisSol currently prints the total simulation time twice via stdout (i.e. to the terminal log), with mostly similar values.
The first value is in the line starting with "Elapsed time (via clock_gettime)",
the second value is printed in the line starting with "Time spent in simulation:", just one line before we see "Simulation done."

The "total time spent in compute kernels" denotes the time in which all CPUs execute some function to advance the computation of the solution.
In particular, it excludes the time spent with MPI communication.

FLOP/s counter
--------------

SeisSol outputs two FLOP and FLOP/s numbers:

* Non-zero FLOP (NZ-FLOP, and NZ-FLOP/s) denote the number of FLOP done towards computing the solution. That is, the number of additions and multiplications while considering the à priori-known sparsity patterns of all involved matrices (thus, e.g. multiplying by a diagonal matrix will be counted as much as a vector product). Note that the degrees of freedom are always assumed to be fully dense.
* Hardware FLOP (HW-FLOP, and HW-FLOP/s) denote the amount of FLOP actually calculated by the CPU. That is, every time we issue an instruction to add or multiply two values, we add a FLOP to the HW-FLOP counter.

Thus, the NZ-FLOP bounded from above by the HW-FLOP. In fact, we could write the HW-FLOP as sum of the NZ-FLOP and the zero-FLOP, i.e. operations where one of the operands is known to be 0 before running the simulation.
This is mostly due to two reasons:

* We do not use the sparsity pattern of a matrix, even though we know it in advance. That can be the case if approximating the sparse matrix by a dense matrix can be implemented more efficiently in hardware by a code generator. Where that may be necessary if usually determined by running an auto tuner (i.e. SeisSol Proxy). By the approximation as a dense matrix, we include some operations where we know one of the operands is to be zero; thus we know the result in advance—therefore these operations are included in the HW-FLOP, but not the NZ-FLOP.
* Padding a matrix, e.g. to match the size of the underlying SIMD registers. The additional values are known to be zero; thus they are only included into the HW-FLOP, not the NZ-FLOP.

From the latter point, it follows that only the HW-FLOP depend on the underlying hardware. To note, both the HW-FLOP and the NZ-FLOP depend on the scenario and the equation system to be solved.

During the simulation (at synchronization points), we only print the HW-FLOP/s. After the simulation has finished, we print both the HW-GFLOP and the NZ-GFLOP, as well as
HW-GLOP/s and NZ-GFLOP/s.

Note that the Dynamic Rupture computation or the Point Sources both are _not_ counted into the HW-/NZ-FLOP numbers at the moment; only the matrix operations do (as used, e.g., during the ADER computation).

Performance
-----------

One important value which we usually publish in our papers is the
performance in "GFLOP/s per node". It comes in the flavor of our two metrics introduced in the previous section:

.. math::

    \text{HW-GFLOP/s per node }= \frac{\text{HW-GFLOP}}{\text{#nodes } \cdot \text{ elapsed-time}}

.. math::

    \text{NZ-GFLOP/s per node }= \frac{\text{NZ-GFLOP}}{\text{#nodes } \cdot \text{ elapsed-time}}

You can compare these values with the publications in order to see if your performance is ok.

Note that, empirically, the "HW-GFLOP/s per node" performance metric is used more often.
