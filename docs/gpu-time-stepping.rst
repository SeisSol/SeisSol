..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _gpu-time-stepping:

Time stepping and halo exchange on GPUs
=======================================

By default, SeisSol advances its time clusters whenever they are ready,
waits for the GPU after each step of a cluster,
and exchanges the halo data between processes with (GPU-aware) MPI.
The options on this page change how the steps get ordered, how the GPU work of the clusters gets enqueued,
and which library exchanges the halo data.
All of them are experimental; the results are meant to stay bitwise identical to the default.

Concepts
~~~~~~~~

Each process runs a set of *clusters*, each taking steps of its own size (local time stepping):

- a *cell cluster* integrates the cells of one time cluster and one layer (interior or copy);
- a *face cluster* computes the dynamic rupture faces of one time cluster and layer;
- a *ghost cluster* stands in for one time cluster of the neighboring processes next to a copy layer,
  and decides when the halo data gets sent and received.

A step consists of a prediction and a correction.
A cluster may take a step once its neighbors have provided the data it needs,
and have read the data it overwrites.
Between two synchronization points (e.g. outputs), the steps are counted in *ticks*,
the steps of the smallest time cluster.
The steps of the largest time cluster form the *super-timesteps*.

Options
~~~~~~~

The options are environment variables.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Meaning
   * - ``SEISSOL_TIMESTEPPING_PLAN=1``
     - Take the steps in a fixed order: by logical time (corrections at the end of their step,
       everything else at its start), corrections before predictions before face work.
       A cluster can then only wait for the halo exchange.
   * - ``SEISSOL_CONCURRENT_CLUSTERS=1``
     - GPU builds only. The clusters only enqueue their GPU work and wait for each other on the GPU (via events),
       instead of waiting for each of their steps. Needs separate scratchpads per layer, which it enables.
   * - ``SEISSOL_SUPERSTEP_GRAPHS=1``
     - Record the GPU work of a super-timestep into a graph, and replay it for all further super-timesteps of the same kind
       (see below). Needs ``SEISSOL_TIMESTEPPING_PLAN=1`` and ``SEISSOL_CONCURRENT_CLUSTERS=1``.
   * - ``SEISSOL_TRANSFER_MODE``
     - How the halo data gets exchanged: ``direct`` (default), ``host``, ``ccl``, ``stream-mpi``, ``shmem``.
       The name ``SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE`` is accepted as well.
   * - ``SEISSOL_EXCHANGE_PER_DIRECTION=1``
     - For ``ccl`` and ``stream-mpi``: one GPU stream per direction between two time clusters,
       instead of one stream for all exchanges in a global order (see below).
       ``SEISSOL_CCL_PER_DIRECTION`` is accepted as well.
   * - ``SEISSOL_SCRATCHPAD_PER_LAYER=1``
     - Give each layer scratchpads of its own, instead of sharing them between all layers.
       The log shows the memory needed either way.
   * - ``SEISSOL_MPI_PERSISTENT=0``
     - Start new MPI requests for each exchange, instead of restarting persistent ones.

Transfer modes
~~~~~~~~~~~~~~

``direct`` and ``host`` exchange the halo data with MPI, from the GPU buffers directly,
or through buffers in host memory for MPI libraries that are not GPU-aware.
The host starts the sends and receives, and tests for their completion.

The remaining modes run the exchange on GPU streams. They need to be enabled at build time:

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Mode
     - Build option
     - Libraries
   * - ``ccl``
     - ``-DCCL=ON``
     - NCCL (CUDA), RCCL (HIP), oneCCL with its C API (oneAPI)
   * - ``stream-mpi``
     - ``-DSTREAM_MPI=MPICH`` or ``-DSTREAM_MPI=CRAY``
     - MPICH 4.1 or newer (``MPIX_Stream``, CUDA and HIP), or the stream-triggered operations of HPE Cray MPICH (``MPIX_Queue``)
   * - ``shmem``
     - ``-DSHMEM=ON``
     - NVSHMEM (CUDA), ROCSHMEM (HIP), Intel SHMEM (oneAPI)

These modes launch the operations of an exchange as one group.
Operations that occupy their stream until the peer has posted the counterpart could wait for each other in a cycle
if every process launched them in its own order.
Hence, by default, all exchanges of a process go through one stream,
ordered by the point in logical time at which their data is complete, which is the same on all processes.
With ``SEISSOL_EXCHANGE_PER_DIRECTION=1``, each direction gets a stream of its own instead;
the groups of different directions may then run at the same time.

``shmem`` puts the data into a staging window of the receiving process in symmetric memory
(as large as the ghost layers of the process with the largest ones),
signals its arrival, and waits until the receiver has cleared the window for the next exchange;
the receiver copies the data from there into its ghost layers.
It always uses the global order.

With ``SEISSOL_CONCURRENT_CLUSTERS=1``, the exchange is ordered on the GPU as well:
the groups wait for the GPU work that produces or last reads their data,
and the clusters wait for the groups they need, all via events.
The host then never waits for the GPU during a super-timestep.

Recording super-timesteps
~~~~~~~~~~~~~~~~~~~~~~~~~

With ``SEISSOL_SUPERSTEP_GRAPHS=1``, the time manager takes the plan one super-timestep at a time.
A super-timestep can be recorded if it is complete up to the synchronization point and no cluster computes on the host.
Two such super-timesteps do the same GPU work if all clusters have the same step sizes.
The first super-timestep of a kind runs as usual, the second one gets recorded into a graph,
and all further ones replay it, while the clusters only keep their books on the host.
All values that change from step to step are read on the GPU: the current time from a clock of each cluster,
and the outputs decide about their samples when their work runs.
Hence, the receivers and fault receivers copy their data to the host in every step in this mode.

Recording needs a GPU that can record graphs, and either no halo exchange (a single process),
or one of the transfer modes on GPU streams in the global order and ``SEISSOL_CONCURRENT_CLUSTERS=1``,
without a communication thread.
If a requirement is missing, a warning says which one, and the super-timesteps run as usual.

Diagnostics
~~~~~~~~~~~

At the end of a run, the log reports the number of halo messages sent and received, summed over all processes;
both have to agree.
With ``SEISSOL_TIMESTEPPING_PLAN=1``, it also reports the super-timesteps:
how many there were, how many ended early at a synchronization point,
how many full ones were free of output samples and of host work, and how many were recorded and replayed.

Status
~~~~~~

The ordering and the decisions of the options above are checked by unit tests,
including simulations of devices and networks that run the enqueued work in any admissible order.
On GPUs, these options have not been validated yet;
compare the results of a run bitwise to the default before relying on them.
The paths for HPE Cray MPICH, ROCSHMEM and Intel SHMEM have not been compiled yet.
