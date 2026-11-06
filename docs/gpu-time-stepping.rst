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
The options on this page change how the steps get ordered, and how the GPU work of the clusters gets enqueued.
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
   * - ``SEISSOL_ACCESS_EVENTS=0``
     - With concurrent clusters, wait for the latest action of each neighbor, instead of only for what the work depends
       on (see below; on by default).
   * - ``SEISSOL_SUPERSTEP_GRAPHS=1``
     - Record the GPU work of a super-timestep into a graph, and replay it for all further super-timesteps of the same kind
       (see below). Needs ``SEISSOL_TIMESTEPPING_PLAN=1`` and ``SEISSOL_CONCURRENT_CLUSTERS=1``.
   * - ``SEISSOL_OUTPUT_THREAD=1``
     - Experimental. Process the outputs of a synchronization point on a thread of their own, while the GPU
       goes on with the next epoch (see below).
   * - ``SEISSOL_HOST_ACTOR=1``
     - Experimental, GPU builds only. Take the samples of the receivers and the fault receivers on the host once
       the GPU is done with what they read, instead of in host functions on streams; also for recorded
       super-timesteps (see below).
   * - ``SEISSOL_OUTPUT_SNAPSHOT=0``
     - With ``SEISSOL_OUTPUT_THREAD=1``: keep no snapshots of the data, and copy the data for the outputs to the host
       before the next epoch, instead of from the snapshots while the GPU goes on (see below).
   * - ``SEISSOL_SUPERSTEP_GRAPH_MODE``
     - Experimental, with ``SEISSOL_SUPERSTEP_GRAPHS=1``. ``replay`` (default): as described below.
       ``always``: besides, record each other super-timestep and launch it once.
       ``rebuild``: record every super-timestep anew and launch it once, replaying none; this measures what the recordings cost.
   * - ``SEISSOL_SUPERSTEP_BRANCHES``
     - With ``SEISSOL_SUPERSTEP_GRAPHS=1``: ``0`` runs the independent branches of the work of a cluster (e.g. the
       neighbor fluxes of the face relations) one after another in the recordings, ``1`` as nodes of their own that can
       run concurrently. With ``0``, the graphs get fewer nodes, which makes them quicker to record and to launch
       (see below). ``1`` by default with CUDA and SYCL, ``0`` with HIP.
   * - ``SEISSOL_TRANSFER_MODE``
     - How the halo data gets exchanged: ``direct`` (default) or ``host``.
       The name ``SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE`` is accepted as well.
   * - ``SEISSOL_SCRATCHPAD_PER_LAYER=1``
     - GPU builds only. Give each layer scratchpads of its own, instead of sharing them between all layers.
       The log shows the memory needed either way.
   * - ``SEISSOL_DEVICE_GRAPH_NODES=0``
     - GPU builds only. Capture all compute graphs from streams, instead of building them node by node where the device allows it.
       The log says which way the graphs get built. No super-timestep gets recorded then.
   * - ``SEISSOL_MPI_PERSISTENT=0``
     - Start new MPI requests for each exchange, instead of restarting persistent ones.
   * - ``SEISSOL_PLAN_INTERIOR_FIRST=0``
     - With ``SEISSOL_TIMESTEPPING_PLAN=1``: correct and do the face work of the copy layers before those of the
       interior layers, as when predicting, instead of after them (on by default: the copy layers wait for the halo,
       and with MPI on the host, the host waits for it in the order of the plan).

Transfer modes
~~~~~~~~~~~~~~

``direct`` and ``host`` exchange the halo data with MPI, from the GPU buffers directly,
or through buffers in host memory for MPI libraries that are not GPU-aware.
The host starts the sends and receives, and tests for their completion.
With concurrent clusters, it starts a send once it has seen the work complete whose data the send reads,
and enqueues the work that needs a receive once the receive has completed.

With concurrent clusters, the work of an action waits for what it depends on: a prediction, which writes the buffers and
derivatives that the neighbors read, for their latest reads (corrections, face work, sends); a correction, which reads
the data of the neighbors, for their latest writes (predictions, face work, receives); the face work for both; the
prediction of a face cluster does nothing and waits for nothing. With ``SEISSOL_ACCESS_EVENTS=0``, every action waits
for the latest action of each neighbor instead, which chains the copy and the interior layer of a cluster one after the
other (their predictions, then their corrections), and all layers of neighboring time clusters as well.

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
and the outputs decide about their samples when their work runs.
Hence, the receivers and fault receivers copy their data to the host in every step in this mode.

The stream runtimes of the clusters and of their outputs build a recording together, node by node:
where the work of one of them waits for that of another, the graph gets an edge, no event.
A replay waits for the work before it, and records one event at its end; the runtimes wait for that event only once they
get work again. So the replays of one super-timestep after another follow each other on one stream, without any wait.

Recording needs a GPU that can build graphs node by node (CUDA 12.3, ROCm 6.3, or SYCL with the oneAPI graph
extension), and no halo exchange, i.e. a single process: the host starts and tests the MPI operations of the
exchange, which a graph cannot replay.
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
On one node with four MI300A (TheRock 23.2.1, November 2026), the recordings with concurrent branches computed other
results than the runs without recordings: on one rank, TPV13 differed in 31 of its 38 outputs, the same way in every run
with two hardware queues (``GPU_MAX_HW_QUEUES``) and another way with four. With one hardware queue, or with the branches one after
another, the outputs were bitwise the same as without recordings; the scratch memory of all branches allocated before
they fork, scratchpads of their own for each layer, the access events, the packet capture of CLR
(``DEBUG_CLR_GRAPH_PACKET_CAPTURE=0``) and an empty node after all leaves of a graph made no difference. The cause is not
known yet; on an RTX 5070, TPV13 on one rank comes out bitwise the same with concurrent branches. So with HIP, the
branches run one after another by default: then TPV13 took 1.75 s on one rank of the MI300A (2.07 s without
recordings), bitwise the same.

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

Work on the host for the GPU
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some work on the host reads what the GPU computes during an epoch: the samples of the receivers
read the unknowns of their cells, and those of the fault receivers the derivatives next to the faults
and the friction state, which the GPU gathers into host memory. By default, a host
function on a stream takes them, which the device runtime calls on a thread of its own once the
gather is done. Such a call is not cheap: on an MI250X, 20–55 µs with ROCm 6.3 and about 500 µs
with ROCm 6.4 and newer (October 2026); and the stream it is on waits for it.

With ``SEISSOL_HOST_ACTOR=1``, the GPU gathers what the work reads into a buffer of host memory of
its own for each step, and an event after the gather stands for it. The thread that enqueues the
work of the GPU checks the events while it goes along (as it does for the halo exchange), and runs
the work of the steps that are done, in their order: by itself, with its OpenMP threads in the loops;
with ``SEISSOL_HOST_TASK_GRAPH=1``, as tasks that the other threads take, in one parallel region per
epoch (see the task graph on the host); or, with ``SEISSOL_OUTPUT_THREAD=1`` (and no task graph), on
the thread of the outputs. The GPU never waits for the host then; there are as many buffers as steps
whose work has not run yet. Before an output reads what the work writes, e.g. the samples at a
synchronization point, all of it runs.

While a super-timestep gets recorded, the GPU gathers that data into device memory of the recording
instead, at addresses the graph keeps. After each launch of the recording, a kernel copies all of it to
a buffer of host memory of its own, and the work goes to the actor, for once the copy is done. So the
recordings contain no host functions, which would hold up the graph while they run. On an RTX 5070 Laptop
GPU, TPV13 with the receivers and fault receivers took 0.63 s instead of 1.84 s with recorded
super-timesteps (0.60 s with tasks), with the same outputs bit for bit.
Checking an event is safe on all backends, as the host reads the data only once the GPU is done with
it, in whatever host memory.
Where the thread that enqueues the work of the GPU is what limits a run, as on LUMI, the work should not run on it:
on two nodes there (MI250X, ROCm 6.3, ``direct``, the time-stepping plan and concurrent clusters), TPV5 took 3.6 s
instead of 3.9 s and TPV13 5.5 s instead of 6.2 s with tasks (``SEISSOL_HOST_TASK_GRAPH=1``), while the work on
the main thread itself gained less and lost in some runs.
On the Radeon 890M (ROCm 6.3), TPV13 with eight receivers took 3.49 s instead of 3.83 s with the
time-stepping plan and concurrent clusters, with the receivers alone on the host actor. On an RTX 5070
Laptop GPU, with the fault receivers as well, 0.67 s instead of 0.76 s (also with tasks, and on the
thread of the outputs); the outputs came out bitwise the same.

Diagnostics
~~~~~~~~~~~

At the end of a run, the log reports the number of halo messages sent and received, summed over all processes;
both have to agree.
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
Scenarios with local time stepping came out bitwise the same as with the default (apart from the energy output,
whose sums vary from run to run anyway), with one process and all options, recorded super-timesteps included,
and with several processes for ``direct`` and ``host``; with SYCL (AdaptiveCpp), with one and two processes and
concurrent clusters. With ``direct``, the time-stepping plan with concurrent clusters and the host actor has run on up
to 128 nodes of Horizon (four B200 each, November 2026) and on LUMI (MI250X, October 2026).
Compare the results of a run bitwise to the default before relying on these options on another system.

On HIP, graphs need Device d2fe2de or newer, which works around two bugs of HIP 7.2:
graph roots that do not wait for the work before the launch, and captures that do not end
when streams waited for each other.
