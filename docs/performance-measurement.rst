..
  SPDX-FileCopyrightText: 2018 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Performance measurement
=======================

Simulation time
---------------

At the end of the simulation, SeisSol prints its time in the lines starting with "Simulation time (total):",
"Simulation time (compute):" and "Simulation time (blocking IO):", each as the mean over the ranks with their range;
the line "Time spent in simulation:", just before "Simulation done.", gives nearly the same total once more.
The compute time is the time spent advancing the clusters, waiting for other ranks included; the blocking IO time is
the time the output took on the main thread.

The "total time spent in compute kernels" denotes the time in which all CPUs execute some function to advance the computation of the solution.
In particular, it excludes the time spent with MPI communication.

.. _run-summary:

Run summary
-----------

What the log says about a run, its performance included, SeisSol also writes for programs to read, next to the
output (``OutputFile`` is the prefix):

* ``<prefix>-summary.json``, the run summary.
* ``<prefix>-epochs.csv``, a row per simulation epoch, i.e. per stretch of simulated time between two synchronization
  points (the times outputs are written at): its simulated time span, its wall time (minimum, mean and maximum over the
  ranks), its compute and blocking IO time (mean and maximum), and the FLOP/s and kernel bytes/s of the epoch,
  summed over the ranks.

Rank 0 writes the summary when the run starts, at every stage of the initialization, every ten seconds or so during
the simulation, and at its end, each time replacing the file at once, so that a reader always finds a complete file.
``status`` is ``completed`` once the run has ended as planned; a run that ended early, by an error or by the time limit
of its job, leaves ``running``, and ``stage`` says how far it got (``init``, ``init.mesh``, ``init.layout``,
``init.model``, ``init.sideConditions``, ``init.io``, ``simulation``; ``finished`` at the end).

Times are in seconds, timesteps and simulated times in simulated seconds. Where a value differs between the ranks, the
summary gives its statistics over them, an object with ``min``, ``max``, ``mean``, ``median``, ``std`` and ``sum``. The
members of the summary are:

.. list-table::
   :header-rows: 1

   * - Member
     - Content
   * - ``schema``
     - the format, ``seissol-run-summary/1``
   * - ``runId``, ``status``, ``stage``, ``started``, ``updated``
     - the ID of the run (see :ref:`env_vars`), how far it got, and when (UTC)
   * - ``build``
     - version, commit, host and device architecture, compiler, build type and the configurations built
   * - ``parallel``
     - ranks, nodes, ranks per node, threads per rank, the MPI library, whether a communication thread runs, the
       ASYNC mode
   * - ``environment``
     - the environment variables of SeisSol, OpenMP, MPI, ASYNC, the device runtimes and the job scheduler
   * - ``invocation``
     - executable, parameter file, working directory, the ``--set`` overrides and the checkpoint the run started from
   * - ``configurations``
     - the configurations the run computes in
   * - ``mesh``
     - cells, copy cells, ghost cells and fault faces per rank
   * - ``lts``
     - the clusters: timestep, update factor, cells and dynamic rupture faces each, and the smallest cell timestep
       before and after the wiggle factor
   * - ``time``
     - ``init``: reading the mesh (``meshRead``), setting it up (``mesh``, which includes reading it) and setting up
       the model (``model``); ``simulation``: the initial output (``initialIo``), the time up to the last epoch
       (``elapsed``), and ``total``, ``compute`` and ``io`` as in the log
   * - ``simulation``
     - start and end time, the time reached, the share of the simulation done, the number of epochs, and whether it
       was aborted (e.g. by an energy criterion)
   * - ``kernels``
     - the time in the compute kernels per rank (``total``), the load imbalance (one minus mean over maximum), and per
       region its time per rank, iterations, samples and the regression of its time over the iterations
   * - ``performance``
     - hardware and non-zero FLOP and kernel bytes, in total, per second and per category

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

During the simulation, we print both HW-FLOP/s and NZ-FLOP/s at synchronization points. After the simulation has finished, we print both the cumulative HW-GFLOP and the NZ-GFLOP, as well as
HW-GLOP/s and NZ-GFLOP/s.

Note that the Dynamic Rupture friction law computation or the Point Sources both are *not* counted into the HW-/NZ-FLOP numbers at the moment; only the matrix operations do (as used, e.g., during the ADER computation).

Kernel Memory Throughput (B/s)
------------------------------

SeisSol outputs one throughput estimation number, namely the kernel memory throughput.

It consists of all non-constant data that is read and written by a kernel;
regardless of where the data resides in memory. As such, it is supposed
to only be an over-estimation over the L2 and main memory bandwidth values.

Constant matrices are assumed to be loaded once at the first kernel launch,
and reside near the compute units (i.e. in L1 or L0, depending on the nomenclature)
during all following executions — thus they are not counted in here.

Due to the structure of SeisSol, the kernel memory throughput will work better for GPUs
than CPUs: on the GPUs (and not on CPUs), a loop over all elements corresponds to a single kernel.
As such, the kernel memory throughput should lie somewhere between L2 and main memory bandwidth
for most kernels.

Performance
-----------

One important value which we usually publish in our papers is the
performance in "TFLOP/s per node" (formerly "GFLOP/s per node"). It comes in the flavor of our two metrics introduced in the previous section:

.. math::

    \text{HW-TFLOP/s per node }= \frac{\text{HW-TFLOP}}{\text{#nodes } \cdot \text{ elapsed-time}}

.. math::

    \text{NZ-TFLOP/s per node }= \frac{\text{NZ-TFLOP}}{\text{#nodes } \cdot \text{ elapsed-time}}

You can compare these values with the publications in order to see if your performance is in line with what to expect.

Note that, empirically, the "HW-TFLOP/s per node" performance metric is used more often; especially in relation to the peak hardware performance.
