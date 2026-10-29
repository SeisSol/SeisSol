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
   * - ``SEISSOL_SUPERSTEP_BRANCHES=0``
     - With ``SEISSOL_SUPERSTEP_GRAPHS=1``: run the independent branches of the work of a cluster (e.g. the neighbor fluxes
       of the face relations) one after another in the recordings, instead of as nodes of their own that can run concurrently.
       The graphs get fewer nodes, which makes them quicker to record and to launch (see below).
   * - ``SEISSOL_TRANSFER_MODE``
     - How the halo data gets exchanged: ``direct`` (default), ``host``, ``ccl``, ``stream-mpi``, ``shmem``, ``actor-mpi``.
       The name ``SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE`` is accepted as well.
   * - ``SEISSOL_EXCHANGE_PER_DIRECTION=1``
     - For ``ccl``, ``stream-mpi`` and ``actor-mpi``: one GPU stream per direction between two time clusters,
       instead of one stream for all exchanges in a global order (see below).
       ``SEISSOL_CCL_PER_DIRECTION`` is accepted as well.
   * - ``SEISSOL_EXCHANGE_LANES``
     - Experimental. For ``ccl``, ``stream-mpi`` and ``actor-mpi`` in the global order: split it into this many lanes
       (1 by default) that run at the same time (see below).
   * - ``SEISSOL_EXCHANGE_SPLIT_BY_NODE=1``
     - Experimental. For ``ccl``, ``stream-mpi``, ``shmem`` and ``actor-mpi``: exchange on the GPU streams only with the
       processes of the own node, and with the others over GPU-aware MPI from the host, as ``direct`` does
       (see below). Not with ``SEISSOL_SUPERSTEP_GRAPHS=1``.
   * - ``SEISSOL_SCRATCHPAD_PER_LAYER=1``
     - GPU builds only. Give each layer scratchpads of its own, instead of sharing them between all layers.
       The log shows the memory needed either way.
   * - ``SEISSOL_DEVICE_GRAPH_NODES=0``
     - GPU builds only. Capture all compute graphs from streams, instead of building them node by node where the device allows it.
       The log says which way the graphs get built. No super-timestep gets recorded then.
   * - ``SEISSOL_MPI_PERSISTENT=0``
     - Start new MPI requests for each exchange, instead of restarting persistent ones.
   * - ``SEISSOL_EXCHANGE_TRANSPORT_EVENTS=1``
     - Experimental. For the modes on GPU streams in the global order with one lane: the work of a cluster waits
       only for the groups of its own transports, instead of all groups launched so far.
   * - ``SEISSOL_ACTOR_MPI_WAIT=kernel``
     - Experimental. For ``actor-mpi``: the streams wait for the host with a kernel of SeisSol that polls with growing
       pauses, instead of the wait of the GPU runtime.
   * - ``SEISSOL_ACTOR_MPI_COPY_ENGINES=1``
     - Experimental. For ``actor-mpi``: copy from and to the host buffers with the copy engines instead of kernels.
   * - ``SEISSOL_ACTOR_MPI_STAGING=0``
     - Experimental. For ``actor-mpi``: MPI sends and receives the regions in GPU memory, instead of host buffers;
       only for MPI libraries whose GPU path does not need the GPU runtime (e.g. HPE Cray MPICH; Open MPI with CUDA hangs).
   * - ``SEISSOL_OWN_QUEUES``
     - Experimental, HIP only, not with ``ccl``, ``stream-mpi`` and ``shmem``. The clusters run their GPU work on
       this many hardware queues of their own (0 by default: none), outside of HIP; all other stream runtimes share one
       more, and the halo exchange of ``actor-mpi`` gets one of its own (with ``actor-mpi``, an odd number gets rounded
       up; see below).
   * - ``SEISSOL_OWN_QUEUES_FACES``
     - Experimental. With ``SEISSOL_OWN_QUEUES``: the face work (dynamic rupture) of the layers on this many own queues
       of its own, next to the work of the cells, instead of on the queues of their cells (0 by default; see below).
   * - ``SEISSOL_OWN_QUEUE_BRANCHES``
     - Experimental. With ``SEISSOL_OWN_QUEUES``: the independent branches of the work of an action (e.g. the face
       relations of the neighboring integration) on this many own queues of their own, shared by all stream runtimes,
       without graphs of single actions (0 by default). A queue of the MI250X runs its kernels one after the other,
       also without barriers in between.
   * - ``SEISSOL_PLAN_INTERIOR_FIRST=0``
     - With ``SEISSOL_TIMESTEPPING_PLAN=1``: correct and do the face work of the copy layers before those of the
       interior layers, as when predicting, instead of after them (on by default: the copy layers wait for the halo,
       and with MPI on the host, the host waits for it in the order of the plan).
   * - ``SEISSOL_OWN_QUEUE_PRIORITIES``
     - Experimental. With ``SEISSOL_OWN_QUEUES``: ``copy`` gives the queues of the copy layers a high priority and those
       of the interior layers a low one (with an even number of queues), ``faces`` the queues of the face work a high
       one; both together as ``copy,faces``.
   * - ``SEISSOL_OWN_QUEUE_LAYOUT=halo``
     - With ``SEISSOL_OWN_QUEUES`` of 2 or more: all copy layers on one queue of a high priority, all interior layers on
       one of a low priority, instead of the copy and interior layer of each time cluster on queues of their own
       (``layers``, the default) as far as there are enough.
   * - ``SEISSOL_ACTOR_MPI_ENGINE=1``
     - Experimental, CUDA only for now. For ``actor-mpi``: a kernel that keeps running during the time stepping does the
       copies of all groups in the order they get ready (see below); ``SEISSOL_ACTOR_MPI_ENGINE_BLOCKS`` (4 by default),
       ``SEISSOL_ACTOR_MPI_ENGINE_FROM_EPOCH`` (the first epoch with it, counted from 0).

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
   * - ``actor-mpi``
     - (none; CUDA and HIP)
     - any MPI (experimental, see below)

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
As the processes launch the groups of different directions in different orders then, operations that wait for their
peer on streams that share a hardware queue can wait for each other in a cycle; NCCL and RCCL require a
consistent order of the operations on different communicators on all processes.
On LUMI (October 2026), RCCL per direction got stuck in the test of the exchange with the four hardware queues per
process that HIP has by default (``GPU_MAX_HW_QUEUES``), and passed it with 16 or 64 of them.

With ``SEISSOL_EXCHANGE_LANES=K``, the global order goes to ``K`` lanes instead: a stream (with a communicator of
its own for ``ccl``) each, of the highest priority, which HIP puts on hardware queues of their own, apart from the
work of the clusters (up to ``GPU_MAX_HW_QUEUES`` of them, minus one for a stream that joins the lanes).
A group goes to a lane by its direction and its exchange, the same on all processes, so each lane has its groups in
the global order on all processes, and the lanes run at the same time.
Ordered on the GPU, the work that depends on the exchange waits for an event that joins all lanes.

With ``SEISSOL_EXCHANGE_SPLIT_BY_NODE=1``, the regions with the processes of other nodes leave the GPU streams:
GPU-aware MPI from the host moves them, once the GPU is done with the work before (as the host sees by polling an
event), while the regions within the node stay on the streams, in the global order. So the latency of the network
does not add up along the order of the exchanges. A send or receive counts as done once MPI has completed its
part, and the streams have launched theirs. For tests, ``SEISSOL_EXCHANGE_RANKS_PER_NODE=N`` takes nodes of ``N``
consecutive ranks each instead of the actual ones.

On two nodes of LUMI (MI250X, ROCm 6.3, eight ranks each, October 2026), TPV5 with the time-stepping plan and
concurrent clusters took 34 s with ``ccl`` in the global order, 4.0 s split by node, and 3.8 s with ``direct``.
Lanes got stuck there: RCCL hung with a second communicator over the network (with the RCCL of TheRock even with
one), while on one node, two and three lanes took 4.3 and 4.8 s (one lane: 4.2 s).

``actor-mpi`` does what ``stream-mpi`` does with any MPI, from the host: a group copies the copy regions of its sends
into buffers in pinned host memory by kernels, and counts that in host memory; a thread of the host sends them once
the count is there, and receives into host buffers as well, counting the completed operations in host memory. The
stream waits for these counts (``cuStreamWaitValue32``, ``hipStreamWaitValue32``) before it copies the received data
into the ghost regions, and before it overwrites the host buffers of a send. So the host never waits for MPI to
enqueue the work of the GPU, and neither the thread nor MPI calls the GPU runtime: the thread that enqueues the work
may block in it meanwhile (e.g. in a launch while the GPU waits), holding locks of the runtime that the copies of a
GPU-aware MPI need as well. Two things on the GPU must not wait for the host then: CUDA must not load a kernel on its
first launch, which may wait for the GPU, so ``actor-mpi`` sets ``CUDA_MODULE_LOADING=EAGER`` unless it is set
already; and a stream that waits blocks the other streams on its hardware queue as well (CUDA has eight by default,
``CUDA_DEVICE_MAX_CONNECTIONS``; on an RTX 5070, eight waiting streams held up all others), which the global order keeps
from waiting for each other in a cycle, as for the other modes.
On a workstation (CUDA, RTX 5070, four ranks on one GPU, October 2026), TPV13 took as long as with ``ccl`` (38 s,
``direct``: 3.1 s; the ranks share the GPU in time slices, and each wait on the GPU costs one).

In the global order on one stream, a group that waits for its data holds up all groups after it. With
``SEISSOL_ACTOR_MPI_ENGINE=1``, the groups only request their copies on the stream instead: a kernel of a few blocks
runs during the time stepping, and does each copy once it is ready, in any order; the work of a cluster then waits only
for the counts of its own transports. At the end of an epoch, it stops once it has done all copies of the epoch.
Before the time stepping, the groups copy on their streams as without it. With CUDA, the first instantiation of some of
the graphs of the clusters waits until the engine has stopped, for a reason not known yet (a test program does not
reproduce it); ``SEISSOL_ACTOR_MPI_ENGINE_FROM_EPOCH=1`` lets the first epoch run without it. On the workstation, an
epoch of TPV13 took 2.5 s with it, against 9.4 s on one stream and 4.4 s on four lanes.
With HIP (ROCm 6.3), the engine hangs: CLR puts all streams, of any priority, and also those with a mask of compute
units, on a few shared hardware queues (``GPU_MAX_HW_QUEUES``), so a kernel that keeps running holds up the streams
that share its queue (13 of 16 streams ran next to it in a test on LUMI). This is also why more lanes do not help there.

On two nodes of LUMI (ROCm 6.3, October 2026; compute time over four or five runs), ``actor-mpi`` took 4.0 to 6.0 s for
TPV5 (``direct``: 3.5 to 4.2 s, ``ccl`` split by node: 3.9 to 5.1 s, ``host``: 3.8 to 4.0 s) and 9.4 to 13.3 s for TPV13
(``direct``: 5.3 to 6.4 s, ``ccl`` split: 6.1 to 7.7 s, ``host``: 5.8 to 6.1 s), bitwise the same as ``direct``. More
lanes, the polling kernel, the copy engines (TPV13: 16 to 18 s), MPI in GPU memory (17 to 20 s), and waiting only for
the own transports did not make it faster; one stream per direction hung. On one node, TPV13 took 12.5 s with ``ccl``
and 12 to 16 s with ``actor-mpi`` (``direct``: 6.2 s): there, the order of the exchange on the GPU costs, not the
network.

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

With concurrent clusters, the work of an action waits for what it depends on: a prediction, which writes the buffers and
derivatives that the neighbors read, for their latest reads (corrections, face work, sends); a correction, which reads
the data of the neighbors, for their latest writes (predictions, face work, receives); the face work for both; the
prediction of a face cluster does nothing and waits for nothing. With ``SEISSOL_ACCESS_EVENTS=0``, every action waits
for the latest action of each neighbor instead, which chains the copy and the interior layer of a cluster one after the
other (their predictions, then their corrections), and all layers of neighboring time clusters as well.

With ``SEISSOL_CONCURRENT_CLUSTERS=1``, the exchange is ordered on the GPU as well:
the groups wait for the GPU work that produces or last reads their data,
and the clusters wait for the groups they need, all via events.
The host then never waits for the GPU during a super-timestep.

Queues of their own (HIP)
~~~~~~~~~~~~~~~~~~~~~~~~~

HIP (its runtime CLR) puts all streams onto a few hardware queues that they share in turn, ``GPU_MAX_HW_QUEUES``
for each priority (4 by default), and creates one more stream for each branch of each graph. A stream that waits
holds up the streams behind it on its queue, and a wait for an event of another queue costs about 9 µs on the GPU
(MI250X, ROCm 6.3; the host needs 12 µs for the launch, the record and the wait); a launch takes 1.7 µs on the host. A
GCD of an MI250X has 24 hardware queues only, and SeisSol has about a hundred streams per process: with 16 queues per
priority, runs hung (the queues that waited kept those they waited for from getting a slot), and from 64 on,
instantiating the graphs ran out of memory. On two nodes of LUMI, two queues per priority were the fastest (TPV5 3.1 to
3.3 s, TPV13 4.9 to 5.2 s, against 3.6 to 3.9 s and 5.5 to 5.9 s with four).

With ``SEISSOL_OWN_QUEUES``, the Device library runs the streams of SeisSol on HSA queues of its own instead
(``device::AbstractAPI::enableOwnQueues()``): the clusters on the given number of queues, all other stream runtimes on
one more, and the halo exchange of ``actor-mpi`` on one of its own; SeisSol leaves HIP one hardware queue per priority
(``GPU_MAX_HW_QUEUES=1``, unless set). The library writes the packets of the queues itself: SeisSol defines
``hipLaunchKernel`` (``device-own-queue-hook``), which hands the kernels launched on such a stream to the library; it
takes the kernel object from the loader and the layout of the arguments from the metadata of the code object. A graph
gets translated into packets when it gets launched on such a stream first, with its temporary allocations moved to
memory of the queue. Events are the marks of a queue (a packet that waits for all work before it, and counts a signal
down); a queue waits for a mark of another one with a packet that the command processor checks (on gfx90a and the
MI300; the RDNA GPUs reject such packets, there a kernel waits). The counters of ``actor-mpi`` are signals, which the
command processor waits for in the same way. A launch then takes 0.7 µs on the host, a wait for another queue about
5 µs, a graph of five small kernels 8 µs (HIP: 91 µs). Events time only on request (``DEVICE_OWN_QUEUE_TIMING=1``,
by a kernel for each event with timing), and host functions run on a thread of the library while their queue waits.
The host never waits to enqueue: packets that find their queue full wait in host memory, until the thread of the
library writes them; with ``actor-mpi``, the queues wait for work that the host enqueues only later.

The fences decide the speed. An acquire of the system in a packet (as HIP has it for its waits) drops the whole L2
cache of the GCD: on the MI250X, a kernel that read 4 MiB again and again took 3.3 times as long next to a queue of such
packets. So a queue waits for another one only with an acquire of the device, which is enough within it; where the
work has to see what the host or MPI wrote into host memory (after a wait for a counter or a host function, and in
copies from host memory), a kernel of the library fences in the shader instead, which drops only the lines that may be
stale (``buffer_invl2``; the hardware keeps those of device memory coherent): 1.1 times as long next to a queue of
those. The marks keep writing back the L2 cache for the system, since MPI reads the device memory past it: with marks
of the device only (``DEVICE_OWN_QUEUE_MARK_SCOPE=agent``), ``direct`` sent stale data.

On two nodes of LUMI (MI250X, ROCm 6.3, October 2026; compute time of two to four runs; the outputs of TPV13 bitwise
the same as with HIP and ``direct``, also with ``actor-mpi`` and with two and four queues), with the time-stepping plan,
concurrent clusters and the host actor (``SEISSOL_HOST_ACTOR=1``):

.. list-table::
   :header-rows: 1
   :widths: 50 25 25

   * - Variant
     - TPV5
     - TPV13
   * - ``direct``, HIP with ``GPU_MAX_HW_QUEUES=2``
     - 3.2 to 3.4 s
     - 4.8 to 4.9 s
   * - ``direct``, one own queue for the clusters
     - 2.9 to 3.1 s
     - 4.5 s
   * - ``direct``, two own queues
     - 2.2 to 2.3 s
     - 3.2 to 3.3 s
   * - ``direct``, four own queues
     - 2.7 s
     - 3.9 s
   * - ``direct``, eight own queues
     - 2.9 s
     - 4.4 to 4.5 s
   * - ``actor-mpi``, HIP with ``GPU_MAX_HW_QUEUES=2`` (without the host actor)
     - 3.9 s
     - 10.6 s
   * - ``actor-mpi``, two own queues
     - 2.2 to 2.3 s
     - 3.3 s
   * - ``actor-mpi``, two own queues, fences of packets for the host (``DEVICE_OWN_QUEUE_HOST_ACQUIRE=packet``)
     - 2.9 s
     - 4.5 s
   * - ``actor-mpi``, two own queues, also waits for other queues of the system (``DEVICE_OWN_QUEUE_WAIT_SCOPE=system``)
     - 3.6 to 3.7 s
     - 5.6 to 5.9 s

With fences of the system for all waits, more own queues made it slower, as with HIP, and one was the fastest (TPV13:
4.8 s, against 5.6 s with two and 7.1 s with eight). With fences of the device, two queues were the fastest, which keep
the copy layers apart from the interior layers (as the ``halo`` layout does, which was as fast with ``actor-mpi``); the
device computed 57 to 59 % of the time for TPV13 and 83 to 84 % for TPV5 then (HIP: 51 % and 73 %; the rank with the
least work 25 % and 29 %). With the waits for what the work depends on (``SEISSOL_ACCESS_EVENTS``, see above), more
queues pay off: TPV13 took 2.7 to 2.9 s on two own queues, 2.4 to 2.5 s on four and 2.4 s on eight (without: 3.9 s on
four), TPV5 1.8 to 1.9 s, 1.6 s and 1.5 to 1.6 s; ``actor-mpi`` on four 2.8 to 3.0 s and 1.6 s (three runs each, the
outputs bitwise the same; with HIP and ``GPU_MAX_HW_QUEUES=2``, they made it 1 to 3 % slower). A trace of the kernels
(``SEISSOL_TRACE_KERNELS``, see :ref:`trace`) showed why: on the ranks with the most work, the queue of the copy layers
ran kernels 85 to 88 % of the time, mostly one right after the other, while kernels of the copy and of the interior
layer that ran at the same time took about as long as alone; the kernels of these small setups leave most of the GPU
idle. On a fixed partition (``SEISSOL_MINISEISSOL=0``: MiniSeisSol weights the partition by a measurement of each rank,
which changes it from run to run), the face work on queues of its own made TPV13 slower, 3.3 s with two of them and 3.0 s
with eight, against 2.5 s on four queues for all (TPV5: 2.1 and 1.9 s against 1.6 s), although none of its loops took
longer: the ranks waited longer for the halo of their neighbors. High priorities for the copy layers or the face work
made no difference, and neither did compute units of their own for the face work (``DEVICE_OWN_QUEUE_CU_MASKS``, e.g.
``-10:0-23;1:24-109`` for queue numbers and ranges of compute units).

Queues hold each other up when they share a pipe of the command processor. On the MI250X, the hardware scheduler gives
the queues of a process four pipes, the i-th queue the pipe i modulo 4 (queues after one that goes move up by one), and
a pipe serves its queues together: a kernel that one queue starts after a long one of another queue on the pipe
completes only after it, and a queue that waits in a packet (for another queue, the host or MPI) holds up the others
on its pipe for about 40 µs at a time (in a test, a queue of short kernels next to four waiting queues lost 40 µs at
every second kernel, next to sixteen at each; kernels that wait instead made it worse, four of them starved a fifth
queue). In SeisSol, the queues of the copy and of the interior layer of different clusters shared pipes: a quarter of
the kernels of the face relations of the copy layers took 50 µs instead of 9. So the library spreads the own queues
over the pipes (``DEVICE_OWN_QUEUE_PIPES=spread``, the default on gfx90a): it finds the pipes with an idle queue on
each, and creates queues until one is on a pipe with the fewest busy own queues; the others stay unused, since the
queues after them would move otherwise. SeisSol declares the queue of the other stream runtimes and that of the
timeline light (``DEVICE_OWN_QUEUE_LIGHT``, queue numbers), which may share a pipe with a busy one;
``DEVICE_OWN_QUEUE_PIPES=check`` only logs the pipes (with ``DEVICE_OWN_QUEUE_VERBOSE=1`` too), ``any`` takes the queues
as they come. On a fixed partition, TPV13 then took 1.9 to 2.0 s on four own queues (2.4 s before), 1.9 s on six and on
eight, 2.2 s on three; TPV5 1.4 s (1.6 s), 1.3 to 1.4 s and 1.6 s; the outputs bitwise the same. The queues of HIP and of
other libraries take pipes as well. Spreading the kernels of a level of a graph over helper queues
(``DEVICE_OWN_QUEUE_GRAPH_LANES=k``, ``DEVICE_OWN_QUEUE_GRAPH_LANE_MIN``) pays off only with pipes left for them: TPV13
2.2 s instead of 2.7 s on two own queues, 2.9 s instead of 1.9 s on four. With
``DEVICE_OWN_QUEUE_GRAPH_BORROW=1``, the lanes go to other own queues on other pipes instead, chosen anew at each
replay among those that are idle, and only while the queue of the graph is idle as well
(``DEVICE_OWN_QUEUE_GRAPH_BORROW_IDLE=0``: whenever the others are); that made no difference on four queues (TPV13
1.94 s, TPV5 1.37 s), and without the condition, which couples the queues, it made them slower (2.03 s and 1.42 s):
which work should run where takes a schedule of the super-steps.

The graphs of the clusters hold most of that latitude: on the rank with the most fault faces of TPV13, the longest path
of dependencies took 87 of 100 ms of a run (with the kernel trace and ``DEVICE_OWN_QUEUE_DAG_LOG``, see
:ref:`trace`), and half of it was the interpolation of the dynamic rupture: per step, eight independent kernels of
about 45 µs, one for each face relation, one after the other on one queue (the four queues computed 58 % of the time).
A simulation of that log over four queues ended the 100 ms in 92 ms with whole graphs placed where they end first,
and in 63 ms with the kernels of the levels of the graphs spread, the longest path first. With
``DEVICE_OWN_QUEUE_SCHEDULE=1``, the graphs spread the kernels of a level over own queues of other pipes where that
ends the level earlier: by the durations each graph measures for its kernels (in a few replays early on, then every
128th), and by the predicted work of the queues (their work as measured or guessed, and that of the queues they wait
for), the longest kernel first onto the queue where it ends first, with a fork and a join of
``DEVICE_OWN_QUEUE_SCHEDULE_FORK_US`` (5 µs) each; a level gets spread only with at least
``DEVICE_OWN_QUEUE_SCHEDULE_MIN_US`` of work (120 µs), and only onto queues alone on their pipe (but for light
ones). On four own queues, TPV13 then took 1.57 to 1.66 s instead of 1.91 to 1.94 s, TPV5 1.06 to 1.09 s instead of
1.36 to 1.38 s, the outputs bitwise the same; with a weight of 1000 for the fault faces in the partition
(``discretization.vertexweightdynamicrupture``, 100 by default; without the schedule, 1.69 s), TPV13 took 1.39 to
1.42 s (TPV5 1.11 s).
Levels of 60 µs of work and less made it as slow as without (forks, joins, and the queues that wait for each other's
work), all levels above 80 µs gave about the same; three own queues took 1.94 and 1.36 s, six 1.67 and 1.18 s, eight
1.88 and 1.36 s (as without). The exchange queue of ``actor-mpi`` waits
nearly all the time and takes a pipe of its own: 2.9 s and 1.6 s on four queues. The clustering of the time steps for
the costs of launches (``LtsClusteringSearch = 'lattice'``, see :ref:`lts_lattice_search`; with the costs of launches and
fills that had made TPV13 14 % faster before) no longer paid off (2.0 s and 1.5 s).

With ``actor-mpi``, the copy layers wait for MPI on the GPU; on queues apart from the interior
layers, these go on meanwhile. On the same queues, TPV13 hung (with one queue for all clusters, before the fences of the
device), so ``actor-mpi`` takes an even number of own queues, which keeps them apart; the priorities of the ``halo``
layout, lanes on queues of their own (``SEISSOL_EXCHANGE_LANES``), transport events and kernel arguments in device
memory (``DEVICE_OWN_QUEUE_KERNARGS=device``) made no difference or made it slower. The exchange engine
(``SEISSOL_ACTOR_MPI_ENGINE=1``) hangs on own queues as well, also on a queue of its own; it is not needed there, since
the command processor waits without a kernel. The library logs its settings once; further options are
``DEVICE_OWN_QUEUE_PACKETS`` (the size of the queues, 131072 packets by default), ``DEVICE_OWN_QUEUE_KERNARG_MIB``
(128) and ``DEVICE_OWN_QUEUE_WAITS=kernel|packet``. The 890M (RDNA 3.5) has 8 hardware queues only; there, runs with
more own queues than fit hung the GPU.

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
Neither can ``actor-mpi`` yet: the host takes up the MPI operations of a group when it enqueues the group.
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
the main thread itself gained less and lost in some runs. Next to RCCL (``ccl``), the OpenMP threads that wait for
tasks should not spin (``OMP_WAIT_POLICY=passive``): they take the cores of the proxy threads of RCCL otherwise
(split by node: 4.9 s instead of 3.9 s).
On the Radeon 890M (ROCm 6.3), TPV13 with eight receivers took 3.49 s instead of 3.83 s with the
time-stepping plan and concurrent clusters, with the receivers alone on the host actor. On an RTX 5070
Laptop GPU, with the fault receivers as well, 0.67 s instead of 0.76 s (also with tasks, and on the
thread of the outputs); the outputs came out bitwise the same.

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
