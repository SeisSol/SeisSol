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

Checkpointing
~~~~~~~~~~~~~

Some environment variables related to checkpointing are described in the :ref:`Checkpointing section <Checkpointing>`.

.. _optimal_environment_variables_on_supermuc_ng:

Optimal environment variables on SuperMUC-NG
--------------------------------------------

On SuperMUC-NG, we recommend using SeisSol with async output in thread mode.
That is SeisSol should be compiled with :code:`-DCOMMTHREAD=ON`, and then run with the environment variables proposed below.
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
