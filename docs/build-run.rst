..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _build_run:

Running SeisSol
===============

At this point, you should have a compiled SeisSol binary ready.
Note that for different scenarios, you will need different binaries, as the underlying equation system varies.

Here, we will only focus on successfully running SeisSol.
And for that, we are going to use the ``SeisSol/Examples`` repository.

Running the proxy
-----------------

To test if your SeisSol installation works, and to get an estimate on its performance,
you can run the SeisSol proxy. To test for full functionality, run

.. code-block:: bash

    SeisSol_proxy_YOUR_CONFIGURATION 1000 10 all

Note that you will need a matching GPU for a SeisSol GPU build.

For benchmarking the performance (make sure you don't do that on a login node), check

.. code-block:: bash

    SeisSol_proxy_YOUR_CONFIGURATION 1000000 50 all

Sometimes, it is useful to only run a specific kernel, e.g. ``ader``.

Running SeisSol on a PC
-----------------------

On your local machine, without having multiple nodes, you may venture into the directory of your choice and run SeisSol.
Note that it may be beneficial in some cases for the performance to set ``OMP_NUM_THREADS`` to the number of cores (not the number of hardware threads) you have.
You should normally not require MPI to run SeisSol with only one process. If you have to, however, make sure to run as ``mpirun -n 1 -bind-to none ./SeisSol`` to still make use of all available cores.
However, if you have two or more GPUs in your PC, you will be required to use N MPI ranks for utilizing all available N GPUs.

Performance considerations
--------------------------

To run SeisSol at full speed, it is highly recommended pinning all threads to CPU cores explicitly, even for GPU builds.

On a cluster, it is advisable to leave one core free from the pinning,
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
You will need to take special care of some masked CPUs, as SeisSol currently does not detect them automatically.
Offline CPUs are detected automatically on Linux systems and avoided in the free CPU masks.

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

IO setup
--------

For the IO, the following environment variables are recommended to be set, when running on a cluster:

.. code-block:: bash

    export XDMFWRITER_ALIGNMENT=8388608
    export XDMFWRITER_BLOCK_SIZE=8388608

    export ASYNC_MODE=THREAD
    export ASYNC_BUFFER_ALIGNMENT=8388608

GPU visibility
--------------

For GPUs, SeisSol should be best launched with one rank per GPU. To select, SeisSol will automatically pick the first visible GPU to it.
However, some systems make all GPUs on a node visible to all processes running on it—potentially resulting in all SeisSol processes
selecting the same GPU. To avoid that, SeisSol provides a launch script given as ``shared/seissol-launch`` which selects the GPU according to the node-local SLURM rank.

Starting SeisSol
----------------

Finally, to run SeisSol, you simply invoke your compiled binary (not the proxy) with a SeisSol parameter file. Like this:

.. code-block:: bash

    ./SeisSol_YOUR_CONFIGURATION parameters.par

If your parameter file is in your launch directory and called ``parameters.par``, you may also leave that parameter away.
