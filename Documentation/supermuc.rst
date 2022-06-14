.. _compile_run_supermuc:


SuperMUC-NG
===========

Setting up GitHub on SuperMuc
-------------------------------

see :ref:`git_behind_firewall`.

Building SeisSol
----------------

1. clone SeisSol including the submodules using 

.. code-block:: bash

  git clone --recursive https://github.com/SeisSol/SeisSol.git
 

2. Load module. Add these lines to .bashrc:

::

  ##### module load for SeisSol
  module load gcc
  module load cmake/3.21.4
  module load python/3.8.11-extended
  module load numactl/2.0.14-intel21
  #To use dependencies preinstalled with spack
  module use /hppfs/work/pr63qo/di73yeq4/myLibs/spack-packages/modules/linux-sles15-skylake_avx512/
  mdule load seissol-env/develop-intel21
  export CC=mpicc 
  export CXX=mpiCC 
  export FC=mpif90

 
3. Alternatively (and for reference), to compile seissol-env on supermuc, follow the procedure below:

.. code-block:: bash

    # load spack
    module load user_spack
    # clone seissol-spack-aid and add the repository
    git clone --branch supermuc_NG https://github.com/SeisSol/seissol-spack-aid.git
    cd seissol-spack-aid
    spack repo add ./spack

    # Manually add a package that cannot be reached by spack
    mkdir -p ~/spack/cache/pkg-config/
    # download package on your local computer and move to supermuc
    wget https://pkgconfig.freedesktop.org/releases/pkg-config-0.29.2.tar.gz
    cp pkg-config-0.29.2.tar.gz ~/spack/cache/pkg-config/

    # install all dependencies of seissol.
    # We specify the intel and intel-mpi version matching preinstalled version on supermuc
    # These can be found with:
    # >spack find intel-mpi
    # >spack compiler list
    spack install seissol-env %intel@21.4.0 ^intel-mpi@2019.12.320
    # or
    spack install seissol-env %gcc@11.2.0 ^intel-mpi@2019.12.320

    # now create a module:
    spack module tcl refresh seissol-env@develop%intel@21.4.0

    #to access the module at start up, add to your ~/.bashrc
    module use $HOME/spack/modules/x86_avx512/linux-sles15-skylake_avx512/
    # change this path to your_custom_path_2_modules if you update ~/.spack/modules.yaml 

Custom install directory for packages and modules can be set with, by changing `` ~/.spack/config.yaml``:

.. code-block:: yaml

    config:
      install_tree: path_2_packages

and `` ~/.spack/modules.yaml``: 

.. code-block:: yaml

    modules:
      default:
        roots:
         tcl: your_custom_path_2_modules

This can be useful to share packages with other user of a supermuc project.

3. Install SeisSol with cmake, e.g. with (more options with ccmake)


.. code-block:: bash

   mkdir build-release && cd build-release
   cmake -DCOMMTHREAD=ON -DNUMA_AWARE_PINNING=ON -DASAGI=ON -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=skx -DPRECISION=double -DORDER=4 -DGEMM_TOOLS_LIST=LIBXSMM,PSpaMM ..
   make -j 48

Note that to use sanitizer (https://en.wikipedia.org/wiki/AddressSanitizer), SeisSol needs to be compiled with gcc.
For that modules and compiler need to be switched:

.. code-block:: bash

    mdule switch seissol-env seissol-env/develop-gcc11

    export CC=mpigcc
    export CXX=mpigxx
    export FC=mpifc

Then cmake on a new build folder.
To enable sanitizer, add ``-DADDRESS_SANITIZER_DEBUG=ON`` to the argument list of cmake, and change the ``CMAKE_BUILD_TYPE`` to ``RelWithDebInfo`` or ``Debug``.

Running SeisSol
---------------

5. Submission file for SeisSol on NG:

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

  # Wall clock limit:
  #SBATCH --time=03:00:00
  #SBATCH --no-requeue

  #Setup of execution environment
  #SBATCH --export=ALL
  #SBATCH --account=<project id>
  #constraints are optional
  #--constraint="scratch&work"
  #SBATCH --partition=general

  #Number of nodes and MPI tasks per node:
  #SBATCH --nodes=40
  #SBATCH --ntasks-per-node=1
  module load slurm_setup
  
  #Run the program:
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
  source /etc/profile.d/modules.sh

  echo 'num_nodes:' $SLURM_JOB_NUM_NODES 'ntasks:' $SLURM_NTASKS 'cpus_per_task:' $SLURM_CPUS_PER_TASK
  ulimit -Ss 2097152
  mpiexec -n $SLURM_NTASKS SeisSol_Release_sskx_4_elastic parameters.par

Accessing PyPI
--------------

Many post-processing scripts of SeisSol require Python dependencies.
We describe how to use pip on SuperMUC at see :ref:`pypi_behind_firewall`.


