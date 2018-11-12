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

| ``export ASYNC_MODE=MPI``
| ``export ASYNC_MPI_COPY=1``
| ``export ASYNC_GROUP_SIZE=64``
| ``export XDMFWRITER_ALIGNMENT=8388608``
| ``export XDMFWRITER_BLOCK_SIZE=8388608``

Note that in the current implementation, at least 2 output nodes have to
be used.

Checkpointing
-------------

-  **SEISSOL_CHECKPOINT_BLOCK_SIZE** Optimize the checkpoints for a
   specific file system block size. Set to 1 to disable the
   optimization. Set to -1 for auto-detection with the SIONlib back-end.
   (default: 1 (MPI-IO, HDF5) or -1 (SIONlib), MPI-IO, HDF5, SIONlib
   back-end only)
-  **SEISSOL_CHECKPOINT_ROMIO_CB_READ** If set, the ``romio_cb_read`` in
   the MPI info object when opening the file. (default: no value, MPI-IO
   and HDF5 back-end only)
-  **SEISSOL_CHECKPOINT_ROMIO_CB_WRITE** See
   *SEISSOL_CHECKPOINT_ROMIO_CB_READ*
-  **SEISSOL_CHECKPOINT_CB_NODES** See
   *SEISSOL_CHECKPOINT_ROMIO_CB_READ*
-  **SEISSOL_CHECKPOINT_ROMIO_DS_READ** See
   *SEISSOL_CHECKPOINT_ROMIO_CB_READ*
-  **SEISSOL_CHECKPOINT_ROMIO_DS_WRITE** See
   *SEISSOL_CHECKPOINT_ROMIO_CB_READ*
-  **SEISSOL_CHECKPOINT_SION_BACKEND** The SIONlib back-end that should
   be used. Should be either *ansi* or *posix*. (default: 'ansi',
   SIONlib back-end only)
-  **SEISSOL_CHECKPOINT_SION_NUM_FILES** Number of SIONlib files.
   (default: 1, SIONlib back-end only)
-  **SEISSOL_CHECKPOINT_SION_COLL_SIZE** The collective size in SIONlib.
   Set to 0 for disabling collective operations. See `SIONlib
   documentation <https://apps.fz-juelich.de/jsc/sionlib/docu/collective_page.html>`__
   for more details. (default: 0, SIONlib back-end only)
-  **SEISSOL_CHECKPOINT_SION_COLL_MODE** The collective mode for
   SIONlib. Should be either *merge* or *normal*. See `SIONlib
   documentation <https://apps.fz-juelich.de/jsc/sionlib/docu/collective_page.html>`__
   for more details. (default: 'merge', SIONlib back-end only)

Optimal environment variables on SuperMuc
-----------------------------------------

Phase 1
~~~~~~~

| Add the following lines prior to invoking SeisSol in your batch file:
| export MP_SINGLE_THREAD=no
| export OMP_NUM_THREADS=16
| export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

Hyperthreading is not efficient on this machine. SeisSol has to be
compiled without commThread ='yes'.

Phase 2
~~~~~~~

.. _for-order-up-to-5-(included):

For order up to 5 (included)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

use the same environement variables as for phase 1, but with 28 threads.

For higher orders
^^^^^^^^^^^^^^^^^

Add the following lines prior to invoking SeisSol in your batch file:

| export MP_SINGLE_THREAD=no
| export OMP_NUM_THREADS=54
| export KMP_AFFINITY=compact,granularity=thread
| export MP_PE_AFFINITY=no

There are 28 CPU on 1 core, but we use hyperthreading. 2 threads are
used for communications. SeisSol has to be compiled with commThread
='yes'.
