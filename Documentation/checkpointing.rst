.. _Checkpointing:

Checkpointing
=============

Checkpointing consists of saving the simulation state at a given time, to allow restarting from that point in case of failure.
It is parametrized within 'output' namelist of the parameter file by setting up the following variables:

.. code:: fortran

   checkPointFile = '../output/check'
   checkPointBackend = 'mpio'
   checkPointInterval = 0.4

| **checkPointFile** defines the path and prefix to the chechpointfile.
| **checkPointBackend** defines the implementation used ('posix', 'hdf5', 'mpio', 'mpio_async', 'sionlib', 'none'). If 'none' is specified, checkpoints are disabled. To use the HDF5, MPI-IO or SIONlib back-ends you need to compile SeisSol with HDF5, MPI or SIONlib respectively.
| **checkPointInterval** defines the (simulated) time interval at which checkpointing is done.  0 (default value) disables checkpointing. When using an asynchronous back-end (mpio_async), you might lose 2 * checkPointInterval of your computation.


If the active checkpoint back-end finds a valid checkpoint during the initialization, it will load it automatically. (You cannot explicitly specify to load a checkpoint)

Hint: Currently only the output of the wave field is designed to work with checkpoints. Other outputs such as receivers and fault output might require an additional post-processing when SeisSol is restarted from a checkpoint.


Checkpointing Environment variables
-----------------------------------

The parallel checkpoint back-ends (HDF5, MPI-IO, SIONlib) support several tuning environment variables:

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


