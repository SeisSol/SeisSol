..
  SPDX-FileCopyrightText: 2019-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _compile_run_supermuc:


SuperMUC-NG
===========

Setting up GitHub on SuperMuc-NG
--------------------------------

see :ref:`pypi_behind_firewall`.

Building SeisSol
----------------

Load modules (including seissol-env). Add these lines to .bashrc:

::

  ##### module load for SeisSol
  module load gcc
  module load cmake/3.21.4
  module load python/3.8.11-extended
  module load numactl/2.0.14-intel21
  #To use dependencies preinstalled with spack
  module use /hppfs/work/pn49ha/ru76tuj2/modules/linux-sles15-skylake_avx512/
  # you need to have access to project pn49ha
  module load seissol-env/develop-intel21-impi-x2b
  export CC=mpicc
  export CXX=mpiCC

Install pspamm (see :ref:`pypi_behind_firewall` for the proxy):

::

    pip3 install git+https://github.com/SeisSol/PSpaMM.git --no-build-isolation --user --proxy http://localhost:DDDDD

(``--no-build-isolation`` is used to circumvent the problem described in https://github.com/SeisSol/PSpaMM/issues/13)

Clone SeisSol including the submodules using

.. code-block:: bash

   git clone --recursive https://github.com/SeisSol/SeisSol.git

Install SeisSol with cmake, e.g. with (more options with ccmake)

.. code-block:: bash

   cd SeisSol
   mkdir build-release && cd build-release
   cmake -DNUMA_AWARE_PINNING=ON -DASAGI=ON -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=skx -DPRECISION=double -DORDER=4 ..
   make -j 48

.. _running_seissol_on_supermuc:

Running SeisSol
---------------

This is an example job submission script for SeisSol on SuperMUC-NG. For your applications, change ``#SBATCH --nodes=``
to the number of nodes you want to run on. A rule of thumb for optimal performance is to distribute your jobs to 1 node per 100k elements. This rule of thumb does not account for potentially shorter queue times, for example when using the test queue or when asking for a large amount of nodes.

::

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
  unset KMP_AFFINITY
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
  srun SeisSol_Release_sskx_4_elastic parameters.par

Accessing PyPI
--------------

Many post-processing scripts of SeisSol require Python dependencies.
We describe how to use pip on SuperMUC at see :ref:`pypi_behind_firewall`.


Using the Sanitizer
-------------------

Note that to use the Sanitizer (https://en.wikipedia.org/wiki/AddressSanitizer), SeisSol needs to be compiled with gcc (or clang but the "static-libasan" argument does not work then).
For that modules and compiler need to be switched:

.. code-block:: bash

    module switch seissol-env seissol-env/develop-gcc11

    export CC=gcc
    export CXX=gxx

Then cmake on a new build folder.
To enable sanitizer, add ``-DADDRESS_SANITIZER_DEBUG=ON`` to the argument list of cmake, and change the ``CMAKE_BUILD_TYPE`` to ``RelWithDebInfo`` or ``Debug``.


Compiling seissol-env spack package
-----------------------------------
For reference, to compile seissol-env on SuperMUC-NG, follow the procedure below:

.. code-block:: bash

    # load spack
    module load user_spack
    # clone seissol-spack-aid and add the repository
    # we use a supermuc specific branch as supermuc spack is too old (0.17.1) for the main branch
    git clone --branch supermuc_NG https://github.com/SeisSol/seissol-spack-aid.git
    cd seissol-spack-aid
    spack repo add ./spack
    # locate externally build pkg-config
    spack external find pkg-config
    # install all dependencies of seissol.
    # We specify the intel and intel-mpi version matching preinstalled version on supermuc-ng
    # These can be found with:
    # >spack find intel-mpi
    # >spack compiler list
    spack install seissol-env %intel ^intel-mpi
    # or
    spack install seissol-env %gcc ^intel-mpi
    # now create a module:
    spack module tcl refresh seissol-env@develop%intel
    #to access the module at start up, add to your ~/.bashrc
    module use $HOME/spack/modules/x86_avx512/linux-sles15-skylake_avx512/
    # change this path to your_custom_path_2_modules if you update ~/.spack/modules.yaml

Custom install directory for packages and modules can be set by changing ``~/.spack/config.yaml``

.. code-block:: yaml

    config:
      install_tree: path_2_packages

and ``~/.spack/modules.yaml``:

.. code-block:: yaml

    modules:
      default:
        roots:
         tcl: your_custom_path_2_modules

This can be useful to share packages with other users of a SuperMUC project.

The seissol-env compilation can also be reduced by adding the python module to ``~/.spack/packages.yaml``:

.. code-block:: yaml

    packages:
      python:
        externals:
        - spec: python@3.8.11
          buildable: False
          modules:
           - python/3.8.11-extended


Compiling the seissol spack package
-----------------------------------

The seissol package is similar to the seissol-env package (it gathers all dependencies of seissol), but also compiles a specific version of seissol itself.
To compile the seissol spack package on SuperMUC-NG, follow the procedure below.

.. code-block:: bash

    # load spack
    module purge
    module load user_spack/23.1.0
    module load intel intel-mkl intel-mpi python/3.10.10-extended

    # install a few python modules (change DDDDD to the value used after RemoteForward in ~/.ssh/config)
    python3.10 -m pip install --upgrade pip --user --proxy http://localhost:DDDDD
    pip3.10 install --upgrade setuptools numpy wheel packaging --user --proxy http://localhost:DDDDD
    pip3.10 install git+https://github.com/SeisSol/PSpaMM.git --no-build-isolation --user --proxy http://localhost:DDDDD

    # clone seissol-spack-aid and add the repository
    # we use a supermuc specific branch as supermuc spack is not fully up to date
    git clone --branch NG https://github.com/SeisSol/seissol-spack-aid.git
    cd seissol-spack-aid
    spack repo add ./spack
    # install a specific version of seissol, and enable python binding enabled for easi
    spack install -j 40 seissol@master convergence_order=4 dr_quad_rule=dunavant equations=elastic precision=single %intel  ^easi +python

    # now create a module:
    spack module tcl refresh $(spack find -d --format "{name}{/hash:5}" seissol)

    #to access the module at start up, add to your ~/.bashrc
    module use $HOME/spack/modules/x86_avx512/linux-sles15-skylake_avx512/
    # change this path to your_custom_path_2_modules if you update ~/.spack/modules.yaml

Custom install directory for packages and modules can be set by changing ``~/.spack/config.yaml``:

.. code-block:: yaml

    config:
      install_tree:
        root: path_2_packages

and ``~/.spack/modules.yaml``:

.. code-block:: yaml

    modules:
      default:
        roots:
         tcl: your_custom_path_2_modules
        tcl:
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
        enable:
        - tcl


This can be useful to share packages with other users of a SuperMUC project.

The number of dependencies to be compiled can be reduced by adding the python module to ``~/.spack/packages.yaml``:

.. code-block:: yaml

    packages:
      python:
        externals:
        - spec: python@3.10.10
          buildable: False
          modules:
           - python/3.10.10-extended
