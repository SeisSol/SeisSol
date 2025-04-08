..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

SuperMUC NG Phase 2
===================

**NOTE: the instructions are still considered highly experimental**

Website: https://doku.lrz.de/supermuc-ng-10745965.html

The nodes consist of (cf. https://doku.lrz.de/hardware-of-supermuc-ng-phase-2-222891050.html ):

- 1× Intel Xeon Platinum 8480+ CPU (Sapphire Rapids), configured with 4 NUMA domains
- 4× Intel Data Center GPU Max 1550 GPUs (Ponte Vecchio), split into two partitions each.

We use ``skx`` as CPU architecture, and ``pvc`` as GPU architecture. As backend (and compiler), we are going to use ``oneapi``.
For compiling SeisSol, we resort to the oneAPI toolkit.

Installing Modules (without Spack)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Firstly, choose a folder you want to install SeisSol to and navigate to it.
Run ``pwd`` and copy the path there. Run the following script there.

.. code-block:: bash

    export SEISSOL_BASE=$(pwd)
    export SEISSOL_PREFIX=$SEISSOL_BASE/local
    export CMAKE_PREFIX_PATH=$SEISSOL_PREFIX:$CMAKE_PREFIX_PATH

    mkdir -p $SEISSOL_PREFIX

Next, we load the necessary modules for our SeisSol build. The default compilers are set automatically.

**WARNING: You may need to install NetCDF by yourself still.**

.. code-block:: bash

    module load python/3.10.12-extended
    module load intel-toolkit/2024.0.0
    module load hdf5/1.14.3-intel24-impi
    module load cmake/3.27.7
    module load ninja/1.11.1

Next, we also start up our Python installation (note that chainforge does not yet support SYCL; hence we resort to gemmforge only. Other tensor compilers are still TBA).

.. code-block:: bash

    python -m venv $SEISSOL_PREFIX
    source $SEISSOL_PREFIX/bin/activate
    pip install setuptools
    pip install numpy
    pip install git+https://github.com/SeisSol/PSpaMM.git
    pip install git+https://github.com/SeisSol/gemmforge.git

Then, we can start installing the modules. For convenience, we also add Ninja as a build tool here first.

.. code-block:: bash

    git clone --branch v1.12.0 --depth 1 https://github.com/ninja-build/ninja.git
    mkdir -p ninja/build
    cd ninja/build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
    make -j 10 install
    cd ../..

The rest of the packages can be installed as usual.

METIS/ParMETIS:

.. code-block:: bash

    wget https://ftp.mcs.anl.gov/pub/pdetools/spack-pkgs/parmetis-4.0.3.tar.gz
    tar -xvf parmetis-4.0.3.tar.gz
    cd parmetis-4.0.3
    sed -i 's/IDXTYPEWIDTH 32/IDXTYPEWIDTH 64/g'  ./metis/include/metis.h
    make config cc=mpicc cxx=mpicxx prefix=$SEISSOL_PREFIX
    make install
    cp build/Linux-x86_64/libmetis/libmetis.a $SEISSOL_PREFIX/lib
    cp metis/include/metis.h $SEISSOL_PREFIX/include
    cd ..

YAML-CPP can be installed as follows:

.. code-block:: bash

    wget https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz
    tar -xf 0.8.0.tar.gz
    mkdir -p yaml-cpp-0.8.0/build
    cd yaml-cpp-0.8.0/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja
    ninja install
    cd ../..

For easi, Eigen and libxsmm, the default instructions suffice.

For ASAGI:

.. code-block:: bash

    git clone --recursive --depth 1 https://github.com/TUM-I5/ASAGI
    mkdir -p ASAGI/build
    cd ASAGI/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja
    ninja install
    cd ../..

For LUA:

.. code-block:: bash

    wget https://www.lua.org/ftp/lua-5.4.6.tar.gz
    tar -xf lua-5.4.6.tar.gz
    cd lua-5.4.6
    make all install INSTALL_TOP=$SEISSOL_PREFIX
    cd ..

For easi (depending on the former two):

.. code-block:: bash

    git clone --recursive --depth 1 https://github.com/seissol/easi
    mkdir -p easi/build
    cd easi/build
    cmake .. -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX -DCMAKE_BUILD_TYPE=Release -GNinja -DASAGI=ON -DLUA=ON -DIMPALAJIT=OFF -DEASICUBE=OFF
    ninja install
    cd ../..

For libxsmm (note that we need 1.17 sharp; the latest main will not work as intended with the generator):

.. code-block:: bash

    git clone --branch 1.17 --depth 1 https://github.com/hfp/libxsmm
    cd libxsmm
    make generator
    cp bin/libxsmm_gemm_generator $SEISSOL_PREFIX/bin
    cd ..

Compiling SeisSol
~~~~~~~~~~~~~~~~~

Finally, it's time to clone SeisSol and build it.

.. code-block:: bash

    git clone --recursive https://github.com/SeisSol/SeisSol.git seissol
    mkdir -p seissol/build
    cd seissol/build
    cmake .. -GNinja -DPRECISION=single -DSYCLCC=dpcpp -DDEVICE_BACKEND=oneapi -DDEVICE_ARCH=pvc -DHOST_ARCH=skx -DORDER=4 -DASAGI=ON -DNUMA_AWARE_PINNING=ON -DCMAKE_INSTALL_PREFIX=$SEISSOL_PREFIX
    ninja

Optionally, you can install SeisSol to ``$SEISSOL_PREFIX``.

Running Jobs
~~~~~~~~~~~~

TBD
