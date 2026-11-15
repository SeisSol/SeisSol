..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _gpu-time-stepping:

Time stepping and halo exchange on GPUs
=======================================

By default, SeisSol advances its time clusters on GPUs whenever they are ready,
lets the clusters only enqueue their GPU work and wait for each other on the GPU (concurrent clusters),
and exchanges the halo data between processes with (GPU-aware) MPI.
The options on this page change how the steps get ordered, and how the GPU work of the clusters gets enqueued.
The results are meant to stay bitwise identical whatever they are set to.

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
   * - ``SEISSOL_CONCURRENT_CLUSTERS=0``
     - GPU builds only. The clusters wait for each of their steps, instead of only enqueueing their GPU work and
       waiting for each other on the GPU (via events; on by default). Concurrent clusters need separate scratchpads
       per layer, which they enable: the log shows the memory that takes.
   * - ``SEISSOL_ACCESS_EVENTS=0``
     - With concurrent clusters, wait for the latest action of each neighbor, instead of only for what the work depends
       on (see below; on by default).
   * - ``SEISSOL_SUPERSTEP_GRAPHS=1``
     - Record the GPU work of a super-timestep into a graph, and replay it for all further super-timesteps of the same kind
       (see below). Needs ``SEISSOL_TIMESTEPPING_PLAN=1`` and concurrent clusters.
   * - ``SEISSOL_OUTPUT_THREAD=1``
     - Experimental. Process the outputs of a synchronization point on a thread of their own, while the GPU
       goes on with the next epoch (see below).
   * - ``SEISSOL_HOST_ACTOR=0``
     - GPU builds only. Take the samples of the receivers and the fault receivers in host functions on streams,
       instead of on the host once the GPU is done with what they read (the host actor, on by default; also for
       recorded super-timesteps, see below).
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
   * - ``SEISSOL_HSA_QUEUES``
     - Experimental, HIP only. The clusters run their GPU work on this many HSA queues, hardware queues outside of HIP
       (0 by default: none), and all other stream runtimes share one more (see below). Only on GPUs where programs
       have run with them (gfx90a, gfx942), unless ``DEVICE_HSA_QUEUE_UNTESTED=1``.
   * - ``SEISSOL_HSA_QUEUES_FACES``
     - Experimental. With ``SEISSOL_HSA_QUEUES``: the face work (dynamic rupture) of the layers on this many HSA queues
       of its own, next to the work of the cells, instead of on the queues of their cells (0 by default).
   * - ``SEISSOL_HSA_QUEUE_BRANCHES``
     - Experimental. With ``SEISSOL_HSA_QUEUES``: the independent branches of the work of an action (e.g. the face
       relations of the neighboring integration) on this many HSA queues of their own, shared by all stream runtimes,
       without graphs of single actions (0 by default).
   * - ``SEISSOL_HSA_QUEUE_PRIORITIES``
     - Experimental. With ``SEISSOL_HSA_QUEUES``: ``copy`` gives the queues of the copy layers a high priority and those
       of the interior layers a low one (with an even number of queues), ``faces`` the queues of the face work a high
       one; both together as ``copy,faces``.
   * - ``SEISSOL_HSA_QUEUE_LAYOUT=halo``
     - Experimental. With ``SEISSOL_HSA_QUEUES`` of 2 or more: all copy layers on one queue of a high priority, all
       interior layers on one of a low priority, instead of the copy and interior layer of each time cluster on queues
       of their own (``layers``, the default) as far as there are enough.

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

A graph that runs only once (``SEISSOL_SUPERSTEP_GRAPH_MODE=always`` or ``rebuild``) tends not to pay off:
instantiating it can take the host longer than enqueuing its work directly.
The more nodes a graph has, the longer the host takes to record it and to launch it;
the log reports the host time of the recordings and of the replays.
Whether concurrent branches (``SEISSOL_SUPERSTEP_BRANCHES=1``) make up for their additional nodes depends on the case.

On MI300A, the recordings with concurrent branches computed other results than the runs without recordings, depending
on the number of hardware queues of HIP (``GPU_MAX_HW_QUEUES``); with one hardware queue, or with the branches one after
another, the outputs were bitwise the same. The cause is not known yet; with CUDA, TPV13 came out bitwise the same with
concurrent branches. So with HIP, the branches run one after another by default.

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
data of the faults), allocated along with it; the log lists their memory per storage.
The host memory of that data is pinned, so that its copies do not block.
With ``SEISSOL_OUTPUT_SNAPSHOT=0``, there are no snapshots, and the data gets copied to the host before the next epoch.
While the copy to the host runs, the kernels that write to the host themselves (e.g. those of the receivers, and the
clocks of the clusters) wait behind it on the link.

With unified memory (``SEISSOL_USM=1``), the snapshots are the host memory of the data that changes:
the GPU computes on memory of its own for that data, and copies it there into the snapshots, which the thread
reads without a further copy; the rest of the data keeps one memory for the host and the GPU.

The outputs came out bitwise the same with the thread and the snapshots as without them.

Work on the host for the GPU
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some work on the host reads what the GPU computes during an epoch: the samples of the receivers
read the unknowns of their cells, and those of the fault receivers the derivatives next to the faults
and the friction state, which the GPU gathers into host memory. With ``SEISSOL_HOST_ACTOR=0``, a host
function on a stream takes them, which the device runtime calls on a thread of its own once the
gather is done. Such a call is not cheap, and the stream it is on waits for it.

By default, the host actor takes them instead: the GPU gathers what the work reads into a buffer of host memory of
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
recordings contain no host functions, which would hold up the graph while they run.
Checking an event is safe on all backends, as the host reads the data only once the GPU is done with
it, in whatever host memory.
Where the thread that enqueues the work of the GPU is what limits a run, the work should not run on it:
tasks (``SEISSOL_HOST_TASK_GRAPH=1``) tend to help there.
The outputs came out bitwise the same with the host actor as without it, also with tasks, on the thread
of the outputs and in recorded super-timesteps.

HSA queues (HIP)
~~~~~~~~~~~~~~~~

HIP (its runtime CLR) puts all streams onto a few hardware queues that they share in turn (``GPU_MAX_HW_QUEUES`` for
each priority, 4 by default), and creates one more stream for each branch of each graph. A stream that waits holds up
the streams behind it on its queue, and a wait for an event of another queue costs the host and the GPU several
microseconds. More hardware queues do not help: a GCD of an MI250X has 24 of them, and SeisSol has about a hundred
streams per process; with 16 hardware queues per priority, runs hung (the queues that waited kept those they waited for
from getting a slot), and with 64, instantiating the graphs ran out of memory.

With ``SEISSOL_HSA_QUEUES``, the Device library runs the streams of SeisSol on HSA queues of its own instead
(``device::AbstractAPI::enableHardwareQueues()``): the clusters on the given number of queues, all other stream runtimes
on one more; SeisSol leaves HIP one hardware queue per priority (``GPU_MAX_HW_QUEUES=1``, unless set). The library
writes the packets of the queues itself: SeisSol defines ``hipLaunchKernel`` (``device-hsa-queue-hook``), which hands
the kernels launched on such a stream to the library; it takes the kernel object from the loader and the layout of the
arguments from the metadata of the code object. A graph gets translated into packets when it gets launched on such a
stream first, with its temporary allocations moved to memory of the queue. Events are the marks of a queue (a packet
that waits for all work before it, and counts a signal down); a queue waits for a mark of another one with a packet
that the command processor checks. Launches and waits for other queues cost the host and the GPU less than with HIP,
graphs far less. Host functions run on a thread of the library while their queue waits, and the host never waits to
enqueue: packets that find their queue full wait in host memory, until the thread of the library writes them.

Events time only on request (``DEVICE_HSA_QUEUE_TIMING=1``, by a kernel for each event with timing). Without, the
loops of the clusters are not timed (see ``SEISSOL_DEVICE_LOOP_STATISTICS`` in :ref:`gpu-env`), and the costs of the
cells cannot be measured for the partition.

On the MI250X and the MI300, the kernel arguments go into device memory that the host writes
(``DEVICE_HSA_QUEUE_KERNARGS=host`` keeps them in host memory, as on the other GPUs), which made the runs faster,
recorded super-timesteps the most, whose kernels read the same arguments again in every replay. A recorded
super-timestep (``SEISSOL_SUPERSTEP_GRAPHS=1``) runs on the one HSA queue it gets launched on, which was slower than
without recordings; ``DEVICE_HSA_QUEUE_GRAPH_ORIGIN=1`` puts its work back on the queues of its streams
(experimental).

The fences decide the speed. An acquire of the system in a packet (as HIP has it for its waits) drops the whole L2
cache of the GCD, which slows down the kernels of all queues. So a queue waits for another one only with an acquire of
the device, which is enough within it; where the work has to see what the host or MPI wrote into host memory (after a
host function, and in copies from host memory), a kernel of the library fences in the shader instead, which drops only
the lines that may be stale (``buffer_invl2``; the hardware keeps those of device memory coherent). The marks keep
writing back the L2 cache for the system, since MPI reads the device memory past it: with marks of the device only
(``DEVICE_HSA_QUEUE_MARK_SCOPE=agent``), ``direct`` sent stale data. With fences of the system for all waits, more HSA
queues made the runs slower, as with HIP; with fences of the device, and with the waits for what the work depends on
(``SEISSOL_ACCESS_EVENTS``, see above), four to eight queues were the fastest on LUMI (MI250X).

Queues hold each other up when they share a pipe of the command processor. On the MI250X, the hardware scheduler gives
the queues of a process four pipes, the i-th queue the pipe i modulo 4 (queues after one that goes move up by one), and
a pipe serves its queues together: a kernel that one queue starts after a long one of another queue on the pipe
completes only after it, and a queue that waits in a packet (for another queue or the host) holds up the others on its
pipe. So the library spreads the HSA queues over the pipes (``DEVICE_HSA_QUEUE_PIPES=spread``, the default on gfx90a):
it finds the pipes with an idle queue on each, and creates queues until one is on a pipe with the fewest busy HSA
queues; the others stay unused, since the queues after them would move otherwise. SeisSol declares the queue of the
other stream runtimes and that of the timeline light (``DEVICE_HSA_QUEUE_LIGHT``, queue numbers), which may share a
pipe with a busy one; ``DEVICE_HSA_QUEUE_PIPES=check`` only logs the pipes, ``any`` takes the queues as they come. The
queues of HIP and of other libraries take pipes as well.

The graphs of the clusters hold most of the latitude: on the rank with the most fault faces of TPV13, the longest path
of dependencies took most of the run (with the kernel trace and ``DEVICE_HSA_QUEUE_DAG_LOG``, see :ref:`trace`), and
half of it was the interpolation of the dynamic rupture: per step, eight independent kernels, one for each face
relation, one after the other on one queue. With ``DEVICE_HSA_QUEUE_SCHEDULE=1``, the graphs spread the kernels of a
level over HSA queues of other pipes where that ends the level earlier: by the durations each graph measures for its
kernels (in a few launches early on, then every 128th), and by the predicted work of the queues, the longest kernel
first onto the queue where it ends first, with a fork and a join of ``DEVICE_HSA_QUEUE_SCHEDULE_FORK_US`` (5 µs) each;
a level gets spread only with at least ``DEVICE_HSA_QUEUE_SCHEDULE_MIN_US`` of work (120 µs), and only onto queues
alone on their pipe (but for light ones). That made TPV5 and TPV13 faster on four HSA queues, and more so with a larger
weight of the fault faces in the partition (*vertexWeightDynamicRupture* in the ``discretization`` section); levels with
less work made it as slow as without (forks, joins, and the queues that wait for each other's work).

The face work on queues of its own (``SEISSOL_HSA_QUEUES_FACES``) made the runs slower (the ranks waited longer for the
halo of their neighbors), and high priorities for the copy layers or the face work made no difference.
The library logs its settings once; further options are ``DEVICE_HSA_QUEUE_PACKETS`` (the size of the queues),
``DEVICE_HSA_QUEUE_KERNARG_MIB`` and ``DEVICE_HSA_QUEUE_WAITS`` (see ``Runtime.cpp`` of the Device library).
The names ``SEISSOL_OWN_QUEUES``, ``SEISSOL_OWN_QUEUE_*`` and ``DEVICE_OWN_QUEUE_*`` of earlier development versions
still work, with a warning.
On other GPUs (``DEVICE_HSA_QUEUE_UNTESTED=1``), check the results first: on the 890M (RDNA 3.5, 8 hardware queues),
runs with more HSA queues than fit hung the GPU.

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
Scenarios with local time stepping came out bitwise the same as with the clusters one after another
(``SEISSOL_CONCURRENT_CLUSTERS=0`` and ``SEISSOL_HOST_ACTOR=0``; apart from the energy output,
whose sums vary from run to run anyway), with one process and all options, recorded super-timesteps included,
and with several processes for ``direct`` and ``host``; with SYCL (AdaptiveCpp), with one and two processes and
concurrent clusters. With ``direct``, the time-stepping plan with concurrent clusters and the host actor has run on
several nodes of Horizon (B200) and of LUMI (MI250X).
On another system, compare the results of a run bitwise to one with these two off before relying on the defaults.
With HSA queues (four, also with the schedule and with the other options; two with the ``halo`` layout), TPV13 and
TPV5 came out bitwise the same as with the streams of HIP on two nodes of LUMI, and with recorded super-timesteps on
one GCD (apart from the energy output).

On HIP, graphs need Device d2fe2de or newer, which works around two bugs of HIP 7.2:
graph roots that do not wait for the work before the launch, and captures that do not end
when streams waited for each other.
