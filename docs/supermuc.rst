..
  SPDX-FileCopyrightText: 2019 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _compile_run_supermuc:


SuperMUC-NG
===========

Discovering Precompiled SeisSol Modules
---------------------------------------

To discover precompiled SeisSol Spack modules on SuperMUC-NG, follow the procedure below.
Update your ``~/.bashrc`` file as follows:

.. code-block:: bash

    # Load Spack
    module purge
    module unload stack
    module load stack/24.5.0
    module unload intel-mpi/2019-intel
    module load intel/2025.1.0
    module load user_spack/24.5.0 python/3.10.12-extended

    # Add the custom SeisSol module path
    module use /hppfs/work/pn49ha/di73yeq4/user_spack24.5/SNG1/modules
    # To access the easi Python module (this points to a symbolic link towards the folder containing the library)
    export PYTHONPATH=/hppfs/work/pn49ha/di73yeq4/user_spack24.5/SNG1/spack_installation/python_wrapper:$PYTHONPATH

After sourcing your ``~/.bashrc``, you should be able to discover the available SeisSol modules, for example:

.. code-block:: bash

    module avail seissol


Setting up GitHub on SuperMuc-NG
--------------------------------

see :ref:`pypi_behind_firewall`.



Installing Python packages with pip
-----------------------------------

The module ``python/3.10.12-extended`` does not support installing packages
using ``pip install --user``. To install additional Python packages, you should
instead create a *virtual environment* with:

.. code-block:: bash

    python -m venv ~/venv-stack24.5-3.10

The virtual environment is then activated with (can be added you ``~/.bashrc``:

.. code-block:: bash

    source ~/venv-stack24.5-3.10/bin/activate

.. _running_seissol_on_supermuc:

Running SeisSol
---------------

This is an example job submission script for SeisSol on SuperMUC-NG. For your applications, change ``#SBATCH --nodes=``
to the number of nodes you want to run on. A rule of thumb for optimal performance is to distribute your jobs to 1 node per 100k elements. This rule of thumb does not account for potentially shorter queue times, for example when using the test queue or when asking for a large amount of nodes.

.. code-block:: bash

  #!/bin/bash
  # Job Name and Files (also --job-name)
  #SBATCH -J <job name>

  #Output and error (also --output, --error):
  #SBATCH -o ./%j.%x.out
  #SBATCH -e ./%j.%x.err

  #Initial working directory:
  #SBATCH --chdir=<work directory>

  #Notification and type
  #SBATCH --mail-type=END
  #SBATCH --mail-user=<your email address>

  #Setup of execution environment
  # note that if you report an issue to LRZ, they will prefer --export=NONE
  #SBATCH --export=ALL
  #SBATCH --account=<project id>
  #SBATCH --no-requeue

  #Number of nodes and MPI tasks per node:
  #SBATCH --partition=general
  #SBATCH --nodes=40
  #SBATCH --time=03:00:00

  #SBATCH --ntasks-per-node=1
  #EAR may impact code performance
  #SBATCH --ear=off

  module load slurm_setup

  #Run the program:
  export MP_SINGLE_THREAD=no
  export OMP_NUM_THREADS=94
  export OMP_PLACES="cores(47)"
  #Prevents errors such as experience in Issue #691
  export I_MPI_SHM_HEAP_VSIZE=32768

  export XDMFWRITER_ALIGNMENT=8388608
  export XDMFWRITER_BLOCK_SIZE=8388608
  export SC_CHECKPOINT_ALIGNMENT=8388608

  export SEISSOL_CHECKPOINT_ALIGNMENT=8388608
  export SEISSOL_CHECKPOINT_DIRECT=1
  export ASYNC_MODE=THREAD
  export ASYNC_BUFFER_ALIGNMENT=8388608

  echo 'num_nodes:' $SLURM_JOB_NUM_NODES 'ntasks:' $SLURM_NTASKS
  ulimit -Ss 2097152
  # NOTE: If you load a SeisSol module, `unset KMP_AFFINITY` must come
  # after loading the module
  unset KMP_AFFINITY

  srun SeisSol_Release_sskx_4_elastic parameters.par

Accessing PyPI
--------------

Many post-processing scripts of SeisSol require Python dependencies.
We describe how to use pip on SuperMUC at see :ref:`pypi_behind_firewall`.

Compiling the seissol Spack package
-----------------------------------

To compile the seissol spack package on SuperMUC-NG, follow the procedure below.

Update your ~/.bashrc to:

.. code-block:: bash

    # load spack
    module purge
    module unload stack
    module load stack/24.5.0
    module unload intel-mpi/2019-intel
    module load intel/2025.1.0
    module load user_spack/24.5.0 python/3.10.12-extended

    # Replace by the port defined in your ~/.ssh/config on your local computer
    if [ "$HOSTNAME" = "login22" ]; then
      export https_proxy=http://localhost:<your_port_for_NG_phase2>
      export http_proxy=http://localhost:<your_port_for_NG_phase2>
    else
      export https_proxy=http://localhost:<your_port_for_NG_phase1>
      export http_proxy=http://localhost:<your_port_for_NG_phase1>
    fi


Clone seissol-spack-aid and add the repository (as spack is not fully up to date on supermuc NG).

.. code-block:: bash

    git clone --branch NG_stack24.5 https://github.com/SeisSol/seissol-spack-aid.git
    cd seissol-spack-aid
    spack repo add ./spack

Adding Python and a few Python modules to the Spack `packages.yaml` configuration file significantly reduces the number of dependencies that need to be compiled.
The `packages.yaml` is discovered with `spack config edit packages`.

.. code-block:: yaml

    packages:
      python:
        buildable: false
        externals:
          - spec: python@3.10.12
            prefix: /lrz/sys/spack/release/sles15.3/24.5.0/views/python/3.10.12-extended/
            modules:
              - python/3.10.12-extended

      py-matplotlib:
        externals:
        - spec: py-matplotlib@3.9.2
          prefix: /lrz/sys/spack/release/sles15.3/24.5.0/views/python/3.10.12-extended/
          buildable: False

      py-scipy:
        externals:
        - spec: py-scipy@1.14.1
          prefix: /lrz/sys/spack/release/sles15.3/24.5.0/views/python/3.10.12-extended/
          buildable: False

      py-numpy:
        externals:
        - spec: py-numpy@2.1.2
          prefix: /lrz/sys/spack/release/sles15.3/24.5.0/views/python/3.10.12-extended/
          buildable: False


tcl modules can be personalized by updating the spack `modules.yaml` config file.
The `modules.yaml` is discovered with `spack config edit modules`.

.. code-block:: yaml

    modules:
      default:
        roots:
         tcl: /hppfs/work/pn49ha/di73yeq4/user_spack24.5/SNG1/modules
         exclude_implicits: false
         auto_load: direct
        tcl:
          hash_length: 0   # <<< This removes the hash suffix
          all:
            suffixes:
              domain_dimension=2: d2
              domain_dimension=3: d3
              polynomial_degree=1: p1
              polynomial_degree=2: p2
              polynomial_degree=3: p3
              polynomial_degree=4: p4
              polynomial_degree=5: p5
              polynomial_degree=6: p6
              convergence_order=3: o3
              convergence_order=4: o4
              convergence_order=5: o5
              convergence_order=6: o6
              equations=elastic: elas
              equations=viscoelastic2: visco
              dr_quad_rule=stroud: stroud
              dr_quad_rule=dunavant: dunav
              precision=single: single
              precision=double: double
              cuda: cuda
              debug: debug
          projections:
            '%oneapi@2021.0:2021.9': '{name}/{version}-oneapi21'
            '%oneapi@2023.0:2023.9': '{name}/{version}-oneapi23'
            '%oneapi@2024.0:2024.9': '{name}/{version}-oneapi24'
            '%oneapi@2025.0:2025.9': '{name}/{version}-oneapi25'
            all: '{name}/{version}'
        enable:
        - tcl


Install a specific version of seissol, and enable python binding enabled for easi:

.. code-block:: bash

    spack install -j 40 seissol@master convergence_order=4 dr_quad_rule=dunavant equations=elastic precision=single %oneapi ^easi +python ^unzip %gcc

Modules can be generated with:

.. code-block:: bash

    spack module tcl refresh $(spack find -d --format "{name}{/hash:5}" seissol)

These modules can be accessed at start up by running `module use`.

.. code-block:: bash

    # change this path to your_custom_path_2_modules if you update ~/.spack/modules.yaml
    module use /hppfs/work/pn49ha/di73yeq4/user_spack24.5/SNG1/modules

Custom install directory for packages and modules can be set by changing the `config.yaml` config file.
The `config.yaml` is discovered with `spack config edit config`.

.. code-block:: yaml

    config:
      install_tree:
        root: /hppfs/work/pn49ha/di73yeq4/user_spack24.5/SNG1/spack_installation

This can be useful to share packages with other users of a SuperMUC project.
