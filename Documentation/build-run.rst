Running SeisSol
===============

At this point, you should have a compiled SeisSol binary ready.
Note that for different scenarios, you will need different.

Here, we will only focus on successfully running SeisSol.
And for that, we are going to use the ``SeisSol/Examples`` repository.

Running the Proxy
-----------------

To test if your SeisSol installation works, and to get an estimate on its performance,
you can run the SeisSol proxy. For that, take

.. code-block::

Running SeisSol on a PC
-----------------------

On your local machine, without having multiple nodes, you may venture into the directory of your choice and run SeisSol.
Note that it may be beneficial in some cases for the performance to set ``OMP_NUM_THREADS`` to the number of cores (not the number of hardware threads) you have.
You should normally not require MPI to run SeisSol with only one process. If you have to, however, make sure to run as ``mpirun -n 1 -bind-to none ./SeisSol`` to still make use of all available cores.
However, if you have two or more GPUs in your PC, you will be required to use N MPI ranks for utilizing all available N GPUs.

Performance Considerations
--------------------------

To run SeisSol properly, it is advisable to do the pinning to the CPU cores explicitly, even for GPU builds.
For now, see the section on environment-variables.rst for a sample configuration. On a cluster, it is advisable to leave one core free from the pinning,
so that it can be used by the MPI communication thread and the IO thread. Generally, the following variables can help:

.. code-block:: bash

    if [ -n "$SLURM_NTASKS_PER_CORE" ]; then
        SEISSOL_HYPERTHREADING = $SLURM_NTASKS_PER_CORE
    else
        SEISSOL_HYPERTHREADING = 1
    fi

    NUM_CORES=$(expr $SLURM_CPUS_PER_TASK / $SEISSOL_HYPERTHREADING)
    NUM_COMPUTE_CORES=$(expr $NUM_CORES - 1)

    export OMP_NUM_THREADS=$NUM_COMPUTE_CORES
    export OMP_PLACES="cores($NUM_COMPUTE_CORES)"
    export PROC_BIND=spread
    export MP_SINGLE_THREAD=no
    unset KMP_AFFINITY

Note that this list is non-exhaustive, and the exact variables are likely cluster-dependent.
You will need to take special care of offline CPUs, as SeisSol currently does not detect them automatically.
Use the environment variable ``SEISSOL_FREE_CPUS_MASK`` to explicitly specify the CPUs that can be used for communication/IO threads.
The variable accepts a comma separated list of elements where an element can be either 1) an integer, or 2) a range of
integers defined as ``[start, end]`` or 3) a comma separated list of integers
surrounded by the curly brackets. The *i*-th list element describes the free cpus
locations for the *i*-th MPI process on the node.

.. code-block:: bash

  # SEISSOL_FREE_CPUS_MASK="(int | range: <int-int> | list: {int,+})+"
  # Examples,
  export SEISSOL_FREE_CPUS_MASK="24,28,32,36"
  export SEISSOL_FREE_CPUS_MASK="24-27,28-31,32-35,36-39"
  export SEISSOL_FREE_CPUS_MASK="{24,25},{28,29},{32,33},{36,37}"

  # Note, it is allowed to mix the formats of list elements. For example,
  export SEISSOL_FREE_CPUS_MASK="24,28-31,{28,29},36"

GPU Visibility
--------------

SeisSol provides a launch script given as ``seissol-launch`` which selects the GPU according to the local SLURM rank.
It is needed due to SeisSol currently only supports one GPU per rank; and the first available device visible is taken automatically.

GPU Memory Management
---------------------



For AMD GPUs, the ``HSA_XNACK`` may need to be set for unified memory to work—if it is supported by your GPU.

Starting SeisSol
----------------


