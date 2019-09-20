Environment Variables
=====================

SeisSol can be tuned with several environment variables:

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

Typical example
^^^^^^^^^^^^^^^

We want to run a 130 nodes job. We dedicate 2 nodes for writing the
asynchronous outputs, the remaining 128 nodes for computing. We then use
a mesh partitioned in 128 regions. In our batch file, we write the
following environment variables:


.. code:: bash

   export ASYNC_MODE=MPI
   export ASYNC_MPI_COPY=1
   export ASYNC_GROUP_SIZE=64
   export XDMFWRITER_ALIGNMENT=8388608
   export XDMFWRITER_BLOCK_SIZE=8388608

Note that in the current implementation, at least 2 output nodes have to
be used.


Checkpointing
~~~~~~~~~~~~~

Some environement variables related to checkpointing are described in the :ref:`Checkpointing section <Checkpointing>`.


Optimal environment variables on SuperMuc
-----------------------------------------

NG
~~

On NG, we recommand using SeisSol with asyn output in thread mode.
That is SeisSol should be compiled with commThread='yes', and then run with the environement variables proposed below.
Also we recommand using hyperthreading capabilities (that is using 96 cpus instead of 48. 2 threads out of 96 are used as communication threads).
Here are some proposed environement variables, to be added prior to invoking SeisSol in your batch file:

.. code:: bash

   #SBATCH --nodes=
   #SBATCH --ntasks-per-node=1 
   #SBATCH --cpus-per-task=96

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


.. _environement_variables_supermuc_phase_2:

Phase 2
~~~~~~~

.. _for-order-up-to-5-(included):

For order up to 5 (included)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| Add the following lines prior to invoking SeisSol in your batch file:

.. code:: bash

   export MP_SINGLE_THREAD=no
   export OMP_NUM_THREADS=28
   export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

SeisSol has to be compiled with commThread ='no'.


For higher orders
^^^^^^^^^^^^^^^^^

Add the following lines prior to invoking SeisSol in your batch file (similarly with NG):

.. code:: bash

   export MP_SINGLE_THREAD=no
   unset KMP_AFFINITY
   export OMP_NUM_THREADS=54
   export OMP_PLACES="cores(27)"

   export XDMFWRITER_ALIGNMENT=8388608
   export XDMFWRITER_BLOCK_SIZE=8388608
   export SC_CHECKPOINT_ALIGNMENT=8388608
   export SEISSOL_CHECKPOINT_ALIGNMENT=8388608
   export SEISSOL_CHECKPOINT_DIRECT=1

   export ASYNC_MODE=THREAD
   export ASYNC_BUFFER_ALIGNMENT=8388608

SeisSol has to be compiled with commThread='yes'.



