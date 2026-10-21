..
  SPDX-FileCopyrightText: 2018 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _env_vars:

Environment variables
=====================

SeisSol can be tuned with several environment variables.

Configuration
-------------

An executable can hold several configurations (see ``EXTRA_CONFIGS`` in the build parameters). ``SEISSOL_CONFIGURATION`` chooses the one a run computes in, by its name, e.g. ``SEISSOL_CONFIGURATION=elastic-linearck-o4-f32-stroud``, for a parameter file that does not name one with ``Configuration`` in its ``equations`` section; a configuration named there takes precedence. Without either, a run takes the first configuration built. So a parameter file written for an executable of a single configuration runs unchanged in an executable of several. The proxy takes the configuration of ``SEISSOL_CONFIGURATION`` as well, unless ``--config`` gives one.

Run ID
------

Every run has an ID, a UUID, which the log prints at the start (``Run ID: ...``) and which every
file the run writes carries: the CSV tables in a ``runId`` column, the HDF5 files as an attribute
``runId`` of their root group, and the text files of the receivers and fault receivers in a line
``# run-id`` of their header. By default, a run takes a new random one. ``SEISSOL_RUN_ID`` gives
one instead, so that whoever starts a run, e.g. a script running a series of them, knows its ID
beforehand; the value has to be a UUID in its canonical form, e.g.
``SEISSOL_RUN_ID=123e4567-e89b-42d3-a456-426614174000``.

Communication thread
--------------------

By default, any SeisSol run with more than one MPI rank will use a communication thread to advance the MPI progress engine.
For that, it is required to leave at least one thread vacant in your OpenMP thread placing map, cf. the SuperMUC-NG example below.

If you do not want to use a communication thread, you may set ``SEISSOL_COMMTHREAD=0`` to make SeisSol poll the progress from time to time.

Task graph on the host
----------------------

By default, every loop over the cells or faces of a cluster is a parallel region of its own, split
among the OpenMP threads evenly; the clusters take their steps one after the other, and all threads
wait at the end of each loop. With ``SEISSOL_HOST_TASK_GRAPH=1``, the main thread decides the steps
of the clusters up to the next synchronization point instead, in one parallel region, and turns the
work of each step into a task, which runs as soon as the work it depends on is done: so the steps of
different clusters run at the same time, and the threads that are done with one loop go on with the
next one. The loops become tasks of ``SEISSOL_HOST_TASK_SIZE`` cells or faces each (32 by default).
The main thread stays at most ``SEISSOL_HOST_ACTION_LIMIT`` steps of the clusters ahead of their work
(256 by default; 0 for no limit): meanwhile, it takes tasks itself. Without a limit, it would decide
all steps up to the next synchronization point first, which costs memory and time with long
intervals between synchronization points.
The halo exchange becomes part of the task graph as well: a send starts once the copy layer is
complete, and the work that reads the ghost layer waits for the data, without the main thread waiting
for either. That pays off the more, the more threads a rank has and the more unevenly they are fast
(e.g. on processors with cores of different kinds); with few threads per rank, it changes little.

The task graph needs OpenMP 5.0 (e.g. GCC 11 or Clang 12 and newer), and it is not available in
builds for a device. With it, a communication thread (see above) starts and tests the exchanges
in the task graph; with ``SEISSOL_COMMTHREAD=0``, the main thread does so between the tasks it takes,
which may be faster where the communication thread shares its core or cache with the workers. That
works with the OpenMP runtimes of LLVM and Intel only: others (e.g. that of GCC) may wait for a task
of an exchange where the main thread creates it, so a thread of its own polls the exchanges there in
any case. The times of the loops in the summary of the run overlap each other in a task graph. The
watchdog (see below) also watches the work done of each cluster then, and at a stall reports the
exchanges of the task graph in flight.

Watchdog
--------

A thread of its own watches the progress of a run, also where the main thread cannot report any
more, e.g. when it waits for the GPU while the halo exchange on the GPU (``ccl``, ``shmem``)
waits for a process that has stopped. If none of the clusters makes progress for
``SEISSOL_WATCHDOG_TIMEOUT`` seconds while the run takes its steps (900 by default), or if anything
else, e.g. the outputs at a synchronization point, takes that long, it reports the phase of the
run, the times of the clusters (also those the GPU has reached), and the clusters report their
state from the threads they act on, the halo clusters with the messages in flight (with MPI, the
rank, the tag and the size of each). Then the run aborts, unless ``SEISSOL_WATCHDOG_ABORT=0``. The
clusters themselves report after the same time without progress. With NCCL and RCCL, the watchdog
also checks for the errors that the library reports asynchronously, e.g. of the network.
``SEISSOL_WATCHDOG=0`` turns it off.

Test of the halo exchange
-------------------------

Before the time stepping, SeisSol tests the halo exchange, through the transfer mode of the run
(see ``SEISSOL_TRANSFER_MODE`` in :ref:`gpu-env`): each copy region gets values that tell the
sending rank, the tag of the region and the position in it, one exchange in each direction moves
them, and each ghost region then has to hold the values of its sender. If a ghost region received
wrong data, each rank concerned names its first such region and position, and the run aborts. A
crash in the library of the exchange, e.g. of an MPI that cannot read device memory, comes right
after the line that announces the test, instead of somewhere in the time stepping. If the exchange
does not complete within ``SEISSOL_WATCHDOG_TIMEOUT`` seconds (and the watchdog is on), the ranks
report the messages still in flight, and the run aborts, unless ``SEISSOL_WATCHDOG_ABORT=0``.
Afterwards, the halo regions are zero again, as before the test. The test costs about one exchange of
the halo, on GPUs plus a few copies of it between the host and the device; ``SEISSOL_HALO_CHECK=0``
turns it off.

Load balancing
--------------

When running with multiple ranks, SeisSol will estimate the node-level performance, enabling better load balancing.
For that, it runs the so-called "Mini SeisSol" benchmark. As its name already hints, it simulates a small test workload on each node;
thus estimating the performance of all nodes relative to each other. The number of elements per node assigned during the partitioning will be resized according to these values.

As a result, the partitioning of runs may become non-deterministic, and the initialization procedure may take a little longer, especially when running only on a single node with multiple ranks.
To disable it, set ``SEISSOL_MINISEISSOL=0``.

In a run with several configurations (see :ref:`several_configurations`), a cell costs according to its configuration. With multiple ranks, Mini SeisSol then also times the time, local and neighbor kernels of each configuration of the run on every rank, and takes the median of the ranks for the cost of a cell of each configuration relative to one of the configuration of the run. The partitioning and the LTS clustering weight each cell with it. Without the measurement (on a single rank, or with ``SEISSOL_MINISEISSOL=0``), the relative cost comes from the hardware FLOPs of the kernels, weighted with the size of a real (4 bytes in single precision, 8 in double precision); that estimate sees nothing of the machine.

Persistent MPI operations
-------------------------

Since SeisSol has a static communication pattern which means, we issue the same MPI transfer requests.
Thus, we use persistent MPI communication which may reduce communication latency.

To disable persistent communication, set ``SEISSOL_MPI_PERSISTENT=0``. Then, SeisSol will use ``MPI_Isend`` and ``MPI_Irecv`` instead. To explicitly enable the persistent communication, set ``SEISSOL_MPI_PERSISTENT=1``. Right now, it is enabled by default.

GPU-Specific Environment Variables
----------------------------------

GPU-specific environment variables are listed on its :ref:`own page <gpu-env>`.

Output
------

* ``SEISSOL_IO_VERTEXFILTER``: whether the coinciding points of an order 0
  output (that is, one written without ``wavefieldvtkorder``, ``surfacevtkorder``
  or ``vtkorder``) are merged into one. On by default; the merging happens
  within a rank, so points on a partition boundary stay duplicated. Also
  accepted under its previous name ``SEISSOL_VERTEXFILTER``.
* ``SEISSOL_IO_META_BLOCK_SIZE``: the size, in bytes, of the blocks the metadata
  of an HDF5 output file is gathered into. Unset, HDF5 keeps its own small
  default. A large value, e.g. ``1048576`` as the previous writer had, keeps the
  many small metadata writes of a file together on a parallel file system, but
  every file then holds at least one such block in full, used or not -- most of
  a file of a few receivers.
* ``SEISSOL_IO_FLUSH``: whether an output file written at every step (a time
  series, the data of an Xdmf output) is flushed after each step. On by default,
  so that the output up to the last step is complete on disk should a run end
  early. A flush syncs the file on every rank, which a parallel file system may
  take its time over with thousands of ranks; turned off, the file is complete
  once the simulation ends.

Two further variables are, for historic reasons, named ``XDMFWRITER``; they are
also accepted with a ``SEISSOL_IO_`` prefix, which is the name to use:

* ``XDMFWRITER_ALIGNMENT``: how far the bulk data of an output is aligned in the
  file. A bulk write straddling a stripe boundary makes more than one storage
  target take part in a single write, which serializes them; what the right
  value is depends on the file system, so there is no default.
* ``XDMFWRITER_BLOCK_SIZE``: no longer read. The previous writer gathered the data
  of an output into blocks of this size before writing them; the MPI-IO hints for
  collective buffering below (``cb_nodes``, ``cb_buffer_size``) take its place.

The MPI-IO hints the output passes through are read the same way, under the
names ROMIO understands and with an ``MPIO_`` in front, so
``SEISSOL_IO_MPIO_ROMIO_DS_WRITE`` sets ``romio_ds_write``. Recognized are
``ind_rd_buffer_size``, ``ind_wr_buffer_size``, ``romio_ds_read``,
``romio_ds_write``, ``cb_buffer_size``, ``cb_nodes``, ``romio_cb_read``,
``romio_cb_write``, ``striping_factor`` and ``striping_unit``. Which of them
help is again a property of the file system.

.. _asynchronous-output:

Asynchronous output
~~~~~~~~~~~~~~~~~~~

In addition to the variables in SeisSol, the
`ASYNC <https://github.com/TUM-I5/ASYNC>`__ library provides some tuning
variables listed in the `wiki <https://github.com/TUM-I5/ASYNC/wiki>`__.

.. _optimal_environment_variables_on_supermuc_ng:

Optimal environment variables on SuperMUC-NG
--------------------------------------------

On SuperMUC-NG, we recommend using SeisSol with async output in thread mode.
Also, we recommend using hyperthreading capabilities (that is, using 96 CPUs instead of 48. 2 threads out of 96 are used as communication threads).
Here are some proposed environment variables to be added before invoking SeisSol in your batch file:

.. code:: bash

   export MP_SINGLE_THREAD=no
   unset KMP_AFFINITY
   export OMP_NUM_THREADS=94
   export OMP_PLACES="cores(47)"

   export XDMFWRITER_ALIGNMENT=8388608
   export XDMFWRITER_BLOCK_SIZE=8388608
   export SC_CHECKPOINT_ALIGNMENT=8388608
   export SEISSOL_CHECKPOINT_ALIGNMENT=8388608
   export SEISSOL_CHECKPOINT_DIRECT=1

   export ASYNC_MODE=THREAD
   export ASYNC_BUFFER_ALIGNMENT=8388608

A complete batch script for SuperMUC-NG can be found in the chapter about :ref:`SuperMUC-NG <running_seissol_on_supermuc>`.

In previous versions of SeisSol, you had to explicitly compile the software with ``-DCOMMTHREAD=ON``. That is not necessary anymore, as
any configuration with more than one MPI rank uses the communication thread by default.
