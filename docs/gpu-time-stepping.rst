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
All clusters count the ticks between the two synchronization times the same way,
rounding up, but ignoring a remainder below :math:`10^{-5}` ticks (the tolerance with which SeisSol compares times);
the last step of each cluster ends exactly at the synchronization point.
The steps of the largest time cluster form the *super-timesteps*.

Options
~~~~~~~

The options are environment variables.
In a build without a GPU, the log warns about the GPU-only options that are set.

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
   * - ``SEISSOL_OUTPUT_THREAD=1``
     - Experimental. Process the outputs of a synchronization point on a thread of their own, while the GPU
       goes on with the next epoch (see below).
   * - ``SEISSOL_OUTPUT_SNAPSHOT=0``
     - With ``SEISSOL_OUTPUT_THREAD=1``: keep no snapshots of the data, and copy the data for the outputs to the host
       before the next epoch, instead of from the snapshots while the GPU goes on (see below).
   * - ``SEISSOL_SUPERSTEP_GRAPH_MODE``
     - Experimental, with ``SEISSOL_SUPERSTEP_GRAPHS=1``. ``replay`` (default): as described below.
       ``always``: besides, record each other super-timestep and launch it once.
       ``rebuild``: record every super-timestep anew and launch it once, replaying none; this measures what the recordings cost.
   * - ``SEISSOL_SUPERSTEP_BRANCHES=0``
     - With ``SEISSOL_SUPERSTEP_GRAPHS=1``: run the independent branches of the work of a cluster (e.g. the neighbor fluxes
       of the face relations) one after another in the recordings, instead of as nodes of their own that can run concurrently.
       The graphs get fewer nodes, which makes them quicker to record and to launch (see below).
   * - ``SEISSOL_TRANSFER_MODE``
     - How the halo data gets exchanged: ``direct`` (default), ``host``, ``ccl``, ``stream-mpi``, ``shmem``.
       The name ``SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE`` is accepted as well.
   * - ``SEISSOL_EXCHANGE_PER_DIRECTION=1``
     - For ``ccl`` and ``stream-mpi``: one GPU stream per direction between two time clusters,
       instead of one stream for all exchanges in a global order (see below).
       ``SEISSOL_CCL_PER_DIRECTION`` is accepted as well.
   * - ``SEISSOL_SCRATCHPAD_PER_LAYER=1``
     - GPU builds only. Give each layer scratchpads of its own, instead of sharing them between all layers.
       The log shows the memory needed either way.
   * - ``SEISSOL_DEVICE_GRAPH_NODES=0``
     - GPU builds only. Capture all compute graphs from streams, instead of building them node by node where the device allows it.
       The log says which way the graphs get built. No super-timestep gets recorded then.
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
With ``ccl``, each direction then also gets a communicator of its own. The first operation towards a peer on a
communicator blocks the process until the peer has joined, to connect them; so all communicators get connected
at the start, one after the other in the same order on all processes, with one exchange of each direction.

``shmem`` puts the data into a staging window of the receiving process in symmetric memory
(as large as the ghost layers of the process with the largest ones),
signals its arrival, and waits until the receiver has cleared the window for the next exchange;
the receiver copies the data from there into its ghost layers.
It always uses the global order.
The signals count the exchanges of each direction. With NVSHMEM and ROCSHMEM, the count of the current exchange is kept
on the GPU: each exchange advances it, and small kernels signal and wait with the device API of the library, reading it
when they run. This needs relocatable device code linked with the device code of the library: with NVSHMEM its static
device library, with ROCSHMEM in a shared library of SeisSol's own (``libseissol-rocshmem-kernels.so``), through which
all calls of ROCSHMEM go; ``-DSHMEM=ON`` sets this up.
Intel SHMEM gets the counts from the host when the operations are enqueued.

With ``SEISSOL_CONCURRENT_CLUSTERS=1``, the exchange is ordered on the GPU as well:
the groups wait for the GPU work that produces or last reads their data,
and the clusters wait for the groups they need, all via events.
The host then never waits for the GPU during a super-timestep.

Recording super-timesteps
~~~~~~~~~~~~~~~~~~~~~~~~~

With ``SEISSOL_SUPERSTEP_GRAPHS=1``, the time manager takes the plan one super-timestep at a time.
A super-timestep can be recorded if it ends before the next synchronization point and no cluster computes on the host.
(The super-timestep that reaches the synchronization point takes the last steps of the clusters, which end exactly there;
their sizes change from one synchronization interval to the next.)
Two such super-timesteps do the same GPU work if all clusters have the same step sizes.
The first super-timestep of a kind runs as usual, the second one gets recorded into a graph,
and all further ones replay it, while the clusters only keep their books on the host.
All values that change from step to step are read on the GPU: the current time from a clock of each cluster,
the exchange counts of ``shmem`` with NVSHMEM or ROCSHMEM,
and the outputs decide about their samples when their work runs.
Hence, the receivers and fault receivers copy their data to the host in every step in this mode.

The stream runtimes of the clusters, of their outputs and of the halo exchange build a recording together, node by node:
where the work of one of them waits for that of another, the graph gets an edge, no event.
A replay waits for the work before it, and records one event at its end; the runtimes wait for that event only once they
get work again. So the replays of one super-timestep after another follow each other on one stream, without any wait.

Recording needs a GPU that can build graphs node by node (CUDA 12.3, ROCm 6.3, or SYCL with the oneAPI graph
extension), and either no halo exchange (a single process),
or ``ccl`` or ``shmem`` with NVSHMEM or ROCSHMEM in the global order and ``SEISSOL_CONCURRENT_CLUSTERS=1``,
without a communication thread.
``stream-mpi`` cannot be recorded: MPICH creates the requests of an enqueued operation on the host
and releases its state in host functions on the stream, which a replay would run again.
Neither can ``shmem`` with Intel SHMEM: a replay would signal and wait for the counts of the recorded exchanges.
If a requirement is missing, a warning says which one, and the super-timesteps (and outputs) run as usual.
Clusters that compute on the host (``SEISSOL_DEVICE_HOST_SWITCH``) keep all super-timesteps from getting recorded;
a warning says so. (Clusters without cells stay on the GPU with any switch point.)

Recording the other super-timesteps as well, each launched once (``SEISSOL_SUPERSTEP_GRAPH_MODE=always``), did not pay off
on a workstation (CUDA, RTX 5070, October 2026): instantiating the graph of a super-timestep took the host 10 to 50 ms
on average, enqueuing its work directly 1 to 2 ms, and a graph that runs once does not make up for that.
With an epoch every 0.05 s, TPV5 and TTPV1 (order 6, single precision) took 0.4 % longer and 0.6 % less time with ``always``
than without recording, and 0.4 % less and 1.0 % more than with the default ``replay``.
Recording every super-timestep anew (``rebuild``) took 20 to 40 % longer.

Capturing a node into a graph costs the more time the more nodes the graph already has, on HIP (ROCm 6.3 to 7.2)
much more than with CUDA: recording a super-timestep of TPV13 (2784 cells, 3 clusters) took 0.04 s on the RTX 5070
and 0.2 s on a Radeon 890M (October 2026). With ROCm 6.3, the host needed 40 ms for each launch of that graph
(2136 kernels), with ROCm 7.2 6 ms; such a graph only pays off where the GPU takes longer for the super-timestep.
The log reports both times for the replays.
With ``SEISSOL_SUPERSTEP_BRANCHES=0``, the graphs of TPV5 (19183 cells) took 9 instead of 51 ms of the host to launch
with ROCm 6.3 on the Radeon 890M, and 1 % less time in all; on the smaller TPV13 mesh, the concurrent branches made up for
their launches, and the run took 3 % (with receivers, 10 %) less time with them.

Outputs while the GPU computes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At a synchronization point, the data that the outputs due there read gets copied to the host,
and only that (e.g. the unknowns for the wave field and the energies; the receivers, the fault
receivers and the ends of epochs without output need none).
With ``SEISSOL_OUTPUT_THREAD=1``, the outputs then take what changes in the next epoch
(e.g. the samples of the receivers, and the data of the elementwise fault output on the GPU), and
process the rest on a thread of their own, while the GPU goes on: the projections, the energies
with their halo exchange and reductions, and the handover to ASYNC.
The thread runs on the CPUs that the OpenMP threads leave free, or else on theirs, which are idle while
the GPU computes, with OpenMP teams of its own, and communicates over a duplicate of the communicator of the run.
The next synchronization point waits for it before it copies to the host again.
An abort the energy output asks for (see the terminator) takes effect once all ranks have processed the same
synchronization points, i.e. an epoch later.
The thread needs all clusters on the GPU, and with unified memory (``SEISSOL_USM=1``) the snapshots (see below);
otherwise, the outputs get processed on the main thread, as without the switch.
A synchronization point where an output cannot be processed on the thread (e.g. of the instantaneous time mirror) gets
processed on the main thread as a whole.

The copy to the host does not hold up the next epoch either: the GPU copies the data within its memory
first, into snapshots, and the thread copies them to the host while the GPU goes on.
For that, the storages keep a snapshot of all the data that changes during the run (e.g. the unknowns and the
data of the faults), allocated along with it; the log lists their memory per storage (in TPV5 on one GPU, 82 MiB).
The host memory of that data is pinned, so that its copies do not block.
With ``SEISSOL_OUTPUT_SNAPSHOT=0``, there are no snapshots, and the data gets copied to the host before the next epoch.
While the copy to the host runs, the kernels that write to the host themselves (e.g. those of the receivers, and the
clocks of the clusters) wait behind it on the link.

With unified memory (``SEISSOL_USM=1``), the snapshots are the host memory of the data that changes:
the GPU computes on memory of its own for that data, and copies it there into the snapshots, which the thread
reads without a further copy; the rest of the data keeps one memory for the host and the GPU.

On a laptop (RTX 5070 Laptop GPU, PCIe link at 7 GB/s, October 2026), TPV5 with an energy output every 0.05 s and the
time-stepping plan with concurrent clusters took 3.19 s without the thread, 2.80 s with it, and 2.69 s with it and the
snapshot (16 % less); with two processes on the GPU, 4.33 s, 4.15 s and 3.83 s (12 % less).
The outputs came out bitwise the same.

Diagnostics
~~~~~~~~~~~

At the end of a run, the log reports the number of halo messages sent and received, summed over all processes;
both have to agree.
With a transfer mode on GPU streams, it also reports how many events of the halo exchange were referenced at once at most,
and how many events their pool has; an event is reused once nobody waits for it any more, so both stay small.
With ``SEISSOL_TIMESTEPPING_PLAN=1``, it also reports the super-timesteps:
how many there were, how many reached a synchronization point (one per synchronization interval and process),
how many of the others (the full ones) were free of output samples and of host work, and how many were recorded and replayed;
with recorded super-timesteps, also the host time they took, by whether they ran directly, were replayed or were recorded,
and how much of the latter went into ending the recordings, which instantiates the graphs.
If a cluster has not taken any action for 15 minutes, the log reports its state and step counts,
and the steps until the synchronization point that it assumes for itself and each of its neighbors.

Status
~~~~~~

The ordering and the decisions of the options above are checked by unit tests,
including simulations of devices and networks that run the enqueued work in any admissible order.
On GPUs, the verification cases pass with all of them.
Beyond that, the options have only been tried on a laptop with an NVIDIA and an AMD GPU (October 2026):
scenarios with local time stepping came out bitwise the same as with the default (apart from the energy output,
whose sums vary from run to run anyway), with one process and all options, recorded super-timesteps included,
and with several processes on one GPU for ``direct``, ``host``, ``ccl`` (NCCL, RCCL) and ``shmem`` (NVSHMEM, ROCSHMEM);
with SYCL (AdaptiveCpp), with one and two processes and concurrent clusters.
None of the transfer modes has run across nodes yet.
The paths for oneCCL and Intel SHMEM have been compiled, but not run;
``stream-mpi`` and the path for HPE Cray MPICH have not been compiled yet.
Compare the results of a run bitwise to the default before relying on these options on another system.

On HIP, graphs need Device d2fe2de or newer, which works around two bugs of HIP 7.2:
graph roots that do not wait for the work before the launch, and captures that do not end
when streams waited for each other.
With CUDA, recording ``ccl`` needs Device 2ec6a42 or newer: CUDA gives all captures into one graph the same id,
which NCCL takes for a single capture, so the NCCL operations get captured into graphs of their own.
SeisSol sets ``NCCL_GRAPH_MIXING_SUPPORT=1`` for ``ccl`` (unless set already), which RCCL needs for recorded super-timesteps.

ROCSHMEM needs version 3.0 or newer (on RDNA GPUs, one whose loads and stores cover their architecture, which those
of 3.2.0 do not). Where ROCSHMEM is built without an MPI at hand (``HAVE_EXTERNAL_MPI`` undefined), it looks up
the predefined objects of Open MPI in the MPI library itself; the version 3.8 of the TheRock nightlies (October 2026)
missed the copies an executable gets by copy relocations there, which are the ones MPI sets up. SeisSol initializes
such a ROCSHMEM from a unique ID shared over MPI instead, so that it bootstraps over TCP.
