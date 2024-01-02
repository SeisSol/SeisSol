Environment Variables
=====================

SeisSol can be tuned with several environment variables.

Communication Thread
--------------------

By default, any SeisSol run with more than one MPI rank will use a communication thread to advance the MPI progress engine.
For that, you will need to leave at least one thread vacant in your OpenMP thread placing map, cf. the SuperMUC-NG example below.

If you do not want to use a communication thread, you may set `SEISSOL_COMMTHREAD=0`; then SeisSol polls on the progress from time to time.

Load Balancing
--------------

When running with multiple ranks, SeisSol will estimate the performance of a node, to enable better load balancing for it.
For that, it runs the so-called "Mini SeisSol" benchmark. As its name already hints at, it simulates a small test workload on each node;
thus estimating the performance of all nodes relative to each other. The number of elements per node assigned during the partitioning will be resized according to these values.

As a result, the partitioning of runs may become non-deterministic, and the initialization procedure may take a little longer; especially when running only on a single node with multiple ranks.
To disable it, set `SEISSOL_MINISEISSOL=0`.

Persistent MPI Operations
-------------------------

Since SeisSol has a static communication pattern (in the sense of: per iteration, we issue the same MPI transfer requests),
we may use persistent MPI communication—it may reduce the communication latency.

You may enable persistent communication by setting `SEISSOL_MPI_PERSISTENT=1`,
and explicitly disable it with `SEISSOL_MPI_PERSISTENT=0`. Right now, it is disabled by default.

Output
------

The wave field and fault output use the
`XdmfWriter <https://github.com/TUM-I5/XdmfWriter>`__. Tuning variables
for the `XdmfWriter <https://github.com/TUM-I5/XdmfWriter>`__ are listed
in the corresponding
`wiki <https://github.com/TUM-I5/XdmfWriter/wiki>`__.

.. _asynchronous-output:

Asynchronous Output
~~~~~~~~~~~~~~~~~~~

In addition to the variables in SeisSol, the
`ASYNC <https://github.com/TUM-I5/ASYNC>`__ library provides some tuning
variables listed in the `wiki <https://github.com/TUM-I5/ASYNC/wiki>`__.

Checkpointing
~~~~~~~~~~~~~

Some environment variables related to checkpointing are described in the :ref:`Checkpointing section <Checkpointing>`.

.. _optimal_environment_variables_on_supermuc_ng:

Optimal environment variables on SuperMUC-NG
--------------------------------------------

On SuperMUC-NG, we recommend using SeisSol with async output in thread mode.
Also, we recommend using hyperthreading capabilities (that is using 96 CPUs instead of 48. 2 threads out of 96 are used as communication threads).
Here are some proposed environment variables, to be added prior to invoking SeisSol in your batch file:

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

In previous versions of SeisSol, you had to explicitly compile the software with `-DCOMMTHREAD=ON`. That is not necessary anymore, as
any configuration with more than one MPI rank uses the communication thread by default.
