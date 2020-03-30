Installation with Spack
=======================

Installation of any HPC application can be tricky. SeisSol is not an exception. 
To considerably alleviate the installation process, we provide few scripts which 
rely on `Spack <https://github.com/spack/spack/wiki>`_, which is a new HPC 
software package manager. 

Spack installs everything from sources and it is the main feature of Spack, 
in contrast to other package managers. In other words, Spack downloads 
source files, compiles them and locates binaries inside of your system. 
You don’t need to have *super-user’s* rights to install packages and 
you can try different versions and flavours of them. However, an installation 
process can take a considerable amount of time due to compilation. Moreover, 
something can go wrong during a compilation process and a user may need to 
read log-files to figure out what exactly went wrong. There are plenty of 
reasons of this. :ref:`Known Issues <spack_known_issues>` is a good starting 
point to have a look in case if you face a problem installing SeisSol 
with Spack.


General Information
-------------------

The installation processes is devided into two parts, 
namely: **seissol-env** and **seissol-utils**.

**seissol-env** - automatically installs all libraries that SeisSol depends on, 
for instance *mpi*, *hdf5*, *netcdf*, *asagi*, etc. 

**seissol-utils** - installs utilities and applications for both pre- and 
post-processing. It is independend of *seissol-env*, however, both can share 
some common dependencies, for example *hdf5*, *netcdf*, etc.


Prerequisites
-------------

First of all, you have to install Spack which you can download from the official 
github `repo <https://github.com/spack/spack.git>`_. Make sure that you have 
*git*, *curl* and *python* installed on your system. Also, be sure that your 
Linux distribution has *build-essential* packages.

.. code-block:: bash

  apt-get install git python curl build-essential


However, most UNIX-kind operating systems come with these packages 
pre-installed. Now, let's install Spack.

.. code-block:: bash

  cd $HOME
  git clone https://github.com/spack/spack.git


Append you *.bashrc* to have Spack avaliable all the time

.. code-block:: bash

  cd spack  
  echo "export SPACK_ROOT=$PWD" >> $HOME/.bashrc
  echo "export PATH=\$SPACK_ROOT/bin:\$PATH" >> $HOME/.bashrc
  cd $HOME


Close and open your terminal to make sure that the changes have been applied. 
For the next step, you need to acquire Spack scripts for SeisSol. 
Download SeisSol from `here <https://github.com/SeisSol/SeisSol>`_ and go to 
the root directory of the application.

.. code-block:: bash

  git clone https://github.com/SeisSol/SeisSol.git
  cd SeisSol


To make SeisSol installation scripts visible inside of Spack, one has 
to add them to Spack repositories. We recommend to install our scripts 
into a separate directory to avoid problems with dangling files inside of 
Spack in case if you decide to delete the current SeisSol repository.


.. code-block:: bash

  cd spack_support
  mkdir build && cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=<install_dir>
  make install
  spack repo add <install_dir>/spack_support


To make sure that everything went well, query avaliable packages in Spack.


.. code-block:: bash

  spack list seissol*
  ==> 2 packages.
  seissol-env  seissol-utils

If you can see an output similar to the one above then we are ready to proceed!

Please, keep in mind that we update installation scripts from time to time. 
Therefore, you have to remove old ones from spack as following:

.. code-block:: bash

  spack repo remove spack_support

Don't forget to add new scripts into the Spack in the same way as we did above.


Getting Started
---------------

One of the main ideas of Spack is to produce a consistent build of your 
software stack, i. e. when everything is compiled with the same compiler suite. 
You may have your preferable compiler suite installed on your system, *intel* 
or *gcc*. If so, you can add them to Spack.

.. code-block:: bash

  spack compiler find <path_to_your_compiler>


However, if you don't have any or you want to try another one, which is not
present in your system, you can install it with Spack. For example, let's 
install *gcc 8.3.0*:

.. code-block:: bash

  spack install gcc@8.3.0


Don't forget to add it to Spack once it is installed:

.. code-block:: bash

  spack compiler find $(spack location -i gcc@8.3.0)


Type the following to see all compilers avaliable for Spack

.. code-block:: bash

  spack compiler list


Environment Modules
-------------------

You can install environment modules to be able to *load* and *unload*
packges, libraries and compilers installed with Spack. 

.. code-block:: bash

  spack bootstrap


After that you can work with the installed software as following:

.. code-block:: bash

  module avail
  module load <package name>
  module list
  module unload <package>
  module purge

You can also look at a list of installed software as following:

.. code-block:: bash

  # the most concise list
  spack find

  # a list of packages with options requested during their instalation
  spack find -v

  # the most detailed list (including install-options of all packages and their deps.)
  spack find -v -d


SeisSol-Env
-----------

The purpose of the script is to install essential packages and libraries for 
SeisSol as well as to install some extra, optional packages that might 
be useful. Here is a list of the essentials:

- hdf5, version=1.8.21
- netcdf-c, version=4.4.0
- libxsmm, version=latest
- pspamm
- memkind, version=latest


*NOTE*: **python3**, **numpy** and **scipy** also belong to the essential 
set and must be on your system to be able to compile SeisSol. However, they 
do not affect run-time performance of SeisSol and most the UNIX-based systems 
have these packages pre-installed. Therefore, installation of these packages 
are optional to save the set-up time. You can trigger an installation of 
them if your system comes without python3 (version=3.5.2), numpy and 
scipy (see, examples). We use the same strategy and reasoning for **cmake** and 
**scons**.

Additionally, a user can customize each individual dependency using 
Spack 
`recursive syntax <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_. 


Options
~~~~~~~

- *asagi* [**default=on**, off] - installs asagi 
- *building_tools* [**default=on**, off] - installs scons and cmake
- *extra_blas* [**default=none**, mkl, openblas, blis] - installs extra blas implementations
- *mpi* [**default=on**, off] - installs an MPI implementation
- *python* [on, **default=off**] - installs python, numpy, scipy and pip

*NOTE*: mpi is a virtual package, a user must specify a concrete implementation
of the standard

Examples
~~~~~~~~

.. code-block:: bash

  # 1. with intel compiler suite
  spack install seissol-env +mpi +asagi %intel@17.0.2 ^intel-mpi@2018.2.199

  # 2. with gcc compiler suite
  spack install seissol-env +mpi +asagi %gcc@8.3.0 ^openmpi@3.1.5

  # 3. with openblas as an extra option
  spack install seissol-env +mpi +asagi extra_blas=openblas %gcc@8.3.0 ^openmpi@3.1.5

  # 4. with a gpu support
  spack install seissol-env +mpi +asagi %gcc@8.3.0 ^openmpi@3.1.5+cuda ^cuda@10.1.243

  # 5. with python, numpy and scipy
  spack install seissol-env +mpi +asagi +python %gcc@8.3.0 ^openmpi@3.1.5


Usage
~~~~~

.. code-block:: bash

  module load seissol-env-develop-<compiler>-<hash>

  # if you compile seissol-env with a compiler installed with Spack
  # you may need to load that compiler as well
  module load <compiler>


After that, you can compile SeisSol using either CMake or 
:ref:`Scons <compiling-seissol>`.


SeisSol-Utils
-------------

By default, the script installs:

- pumgen (without a Simmetrix support)
- gmsh (without a GPU support)
- gmsh2gambit
- cube_c
- rconv
- SeisSol Cookbook, which contains some examples to run

As in case of *seissol-env*, you need **scons** and, therefore, **python3** for 
compiling. However, installation of these packages is optional to save 
the set-up time.


Options
~~~~~~~

- *benchmarks* [on, **default=off**] - installs SeisSol benchmarks. Make sure that you have access to the SeisSol LRZ-gitlab account.
- *building_tools* [on, **default=off**] - installs scons and as a result python and pip
- *gmsh_gui* [on, **default=off**] - enables gui support for gmsh
- *paraview* [on, **default=off**] - installs Paraview for visualization

Examples
~~~~~~~~

.. code-block:: bash

  # 1. essential packages compiled with gcc compiler suite
  spack install seissol-utils %gcc@8.3.0

  # 2. with benchmarks and gmsh gmsh GUI
  spack install seissol-utils+gmsh_gui+benchmarks %gcc@8.3.0

  # 3. with gmsh GUI, paraview and scons
  spack install seissol-utils+gmsh_gui+paraview+building_tools %gcc@8.3.0

  # 4. essential packages with simmetrix support for pumgen
  spack install seissol-utils %gcc@8.3.0 ^pumgen+simmetrix_support 

Usage
~~~~~

.. code-block:: bash

  module load seissol-utils-develop-<compiler>-<hash>

  # to access the Cookbook
  cd $COOKBOOK

  # to access the Benchmakrs
  cd $BENCHMAKRS



Tips and Tricks
---------------

.. _spack_known_issues:

1. Spack builds the entire dependency graph before compiling and installing. 
The graph includes all libs and packages which are necessary to build your 
application, including packages like: *tar, gzip, zlib,  autoconf, 
cmake, automake, pkgconf, m4, ncurses, etc*. Packages like these do not 
affect performance of your application but help Spack to install it. 
Therefore, it is not necessary to install them again and again. You can 
install such  packages only once and mark them as Default 
`(External) <https://spack-tutorial.readthedocs.io/en/latest/tutorial_configuration.html#external-packages>`_.
and Non-Buildable. It can speed-up installation of SeisSol-Env and SeisSol-Utils 
considerably. You will need to modify and edit **~/.spack/packages.yaml** file.


Known Issues
------------

1. Spack is a really live project with dozens of commits per day. It is 
difficult for us to keep the same pace with Spack. A new version of Spack
may not work because of new added features what we may not be aware of. 
Therefore, it may be necessary to use an older version of Spack. You
can simply do it by moving the HEAD of your locally installed Spack
repository to an old commit:

.. code-block:: bash

  cd $SPACK_ROOT
  git checkout <a previous SPACK commit>



2. You may need to reload **setup-env.sh** script if you cannot see 
packages in the module system right after their installation.

.. code-block:: bash

  source $SPACK_ROOT/share/spack/setup-env.sh


3. Some low-level packages are sensitive to your environment variables and 
small syntactic mistakes can lead to weird compilation errors. Please, check 
your environment variables in advance to avoid it. Make sure that you don't 
have trailing or leading **colons and dots** in PATH, LD_LIBRARY_PATH, 
C_INCLUDE_PATH, etc.


4. Some compilers, especially new ones, are not always able to successfully 
install all SeisSol software stack. If it is a case you can try the 
installation process again using an older version of your compiler.


5. Spack is an HPC package manager. Most of HPC systems have a fast-access 
file storage attached to **/tmp** directory to handle temporary files as 
fast as possible. Spack knows about it and takes advantage out of it. 
By default, Spack use **/tmp** for compiling, building and caching your 
binaries.  If you software stack is relatively huge and you would like 
to have multiple versions of your software stack compiled with different 
‘flavours’ this directory can quickly exhaust the memory space allocated 
for your system. Usually, your home directory is attached to a slower but 
bigger storage-drive and sometimes it is better to change the default 
Spack behavior. You will have to modify **~/.spack/config.yaml** file. 
For example:

.. code-block:: bash

  cat ~/.spack/config.yaml
  config:                                                                                                               
      build_stage:                                                                                                      
          - ~/.tmp_build                                                                                                
          - ~/.spack/stage
