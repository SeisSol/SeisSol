.. _compile_run_supermuc:


SuperMUC-NG
===========

Accessing github
----------------

SuperMUC restricts access to outside sources and thus does not allow connections to https servers. 
Nevertheless, GitHub can be used if remote port forwarding is correctly set.
Here, we described the procedure to set up such port forwarding.


1. On your local machine, add to your ~/.ssh/config the following lines:

::

  Host supermucNG
     Hostname skx.supermuc.lrz.de
     User <Your Login>    
     RemoteForward ddddd github.com:22

where ddddd is an arbitrary 5-digital port number, smaller than 65535.
  
2. Use the following command to login onto SuperMUC-NG:

::

  ssh supermucNG 
  
Add the following lines to your ~/.ssh/config (on supermucNG):

:: 

  Host github.com
     HostName localhost
     User git
     Port ddddd
    
With ddddd the same port number as before.

3. Create SSH key by typing (use a non-empty passphrase, not too long as you will need to type it often)

::

  ssh-keygen -t rsa 

4. Go to https://github.com/settings/ssh, add a new SSH key, and paste the public SSH key you just created (the content of ~/.ssh/id_rsa.pub on supermucNG). You should now be able to clone SeisSol including the submodules using:

::

  git clone git@github.com:SeisSol/SeisSol.git

Pay attention to the change in the git address ('https://github.com/' is now replaced by 'git@github.com:'). 
If it works, you will see several lines, for example: 

::

  Cloning into 'SeisSol'...
  remote: Enumerating objects: 25806, done.
  remote: Counting objects: 100% (4435/4435), done.
  remote: Compressing objects: 100% (1820/1820), done.
  remote: Total 25806 (delta 2972), reused 3710 (delta 2551), pack-reused 21371
  Receiving objects: 100% (25806/25806), 110.50 MiB | 9.79 MiB/s, done.
  Resolving deltas: 100% (19382/19382), done.


Building SeisSol
----------------

1. clone SeisSol including the submodules using 

::

  git clone git@github.com:SeisSol/SeisSol.git
  cd SeisSol
  git submodule update --init
 

2. Load module. Could add these lines to .bashrc (changing the order and adding additional modules may prevent a successful compilation):

::

  ##### module load for SeisSol
  module load gcc
  module load cmake/3.21.4
  module load python/3.8.11-extended
  module load libszip/2.1.1
  module load netcdf-hdf5-all/4.7_hdf5-1.10-intel21-impi
  module load numactl/2.0.14-intel21
  module load yaml-cpp/0.7.0-intel21


3. Install eigen, metis, parmetis, libxsmm, PSpaMM, easi and ASAGI or rely on our precompiled libraries

We installed all these libraries in /hppfs/work/pr63qo/di73yeq4/myLibs/SeisSol_dependencies_intel
Simply add the following line to your ~/.bashrc file:

.. code-block:: bash

    export SeisSolDepFolder=/hppfs/work/pr63qo/di73yeq4/myLibs/SeisSol_dependencies_intel
    export PATH=$SeisSolDepFolder/bin/:$PATH
    export PKG_CONFIG_PATH=$SeisSolDepFolder/lib/pkgconfig:$PKG_CONFIG_PATH
    export LD_LIBRARY_PATH=$SeisSolDepFolder/lib:$LD_LIBRARY_PATH

Alternatively, to install by yourself, following the instructions below, you need to change the first line (SeisSolDepFolder) to your home directory (because that where libraries are installed when following the documentation) :

.. code-block:: bash

    export SeisSolDepFolder=~

Then, follow :ref:`installing_eigen3`, :ref:`installing_parmetis`, :ref:`installing_libxsmm`, :ref:`installing_pspamm`, :ref:`installing_ASAGI` and `Installing easi <https://easyinit.readthedocs.io/en/latest/getting_started.html>`_.
Note that ASAGI needs to be compiled before easi.

4. Install SeisSol with cmake, e.g. with (more options with ccmake)

.. code-block:: bash

   mkdir build-release && cd build-release
   CC=mpicc CXX=mpiCC FC=mpif90  cmake -DCMAKE_PREFIX_PATH=$SeisSolDepFolder -DCOMMTHREAD=ON -DNUMA_AWARE_PINNING=ON -DASAGI=ON -DCMAKE_BUILD_TYPE=Release -DHOST_ARCH=skx -DPRECISION=double -DORDER=4 -DCMAKE_INSTALL_PREFIX=$(pwd)/build-release -DGEMM_TOOLS_LIST=LIBXSMM,PSpaMM -DPSpaMM_PROGRAM=$SeisSolDepFolder/bin/pspamm.py ..
   make -j 48

Note that to use sanitizer (https://en.wikipedia.org/wiki/AddressSanitizer), SeisSol needs to be compiled with gcc.
For that modules and compiler need to be switched:

::

    module load netcdf-hdf5-all/4.7_hdf5-1.10-gcc11-impi
    module load numactl/2.0.14-gcc11
    module load yaml-cpp/0.7.0

    export CC=mpigcc
    export CXX=mpigxx
    export FC=mpifc

Then cmake (without ``CC=mpicc CXX=mpiCC FC=mpif90``) on a new build folder.
easi (and all its dependencies) also needs to be build with gcc compilers.
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
We describe how to use pip on SuperMUC in the following.


1. On your local machine in ~/.ssh/config add the following `RemoteForward` line:

::

    Host supermucNG
        ...
        RemoteForward ddddd localhost:8899

where ddddd is an arbitrary port number with 5 digits.
(This number should be different from port number used in other RemoteForward entries.)

2. Install proxy.py on your local machine.

::

    pip install --upgrade --user proxy.py

3. Start proxy.py on your local machine. (And keep it running.)


::

    ~/.local/bin/proxy --port 8899

4. Login to SuperMUC with `ssh supermucNG`. Pip can be used with

::

    pip install <package name> --user --proxy localhost:ddddd

where ddddd is your arbitrary port number.
