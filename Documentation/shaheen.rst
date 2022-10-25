.. _compile_run_shaheen:


Shaheen
=======


Add these lines to ``~/.bashrc``:

::

    ##### module load for SeisSol
    module unload PrgEnv-cray
    module load PrgEnv-gnu
    module load python
    module load cmake
    module unload intel
    # to use preinstall seissol-env module
    module use /project/k1587/ulrich/spack/share/spack/modules/cray-cnl7-ivybridge/
    module load seissol-env-develop-gcc-11.2.0-rckyrcj
    export CC=cc
    export CXX=CC
    export FC=ftn

This will load a preinstalled seissol-env module.

Alternatively (and for reference), to compile seissol-env on shaheen, follow the procedure below:

.. code-block:: bash

    git clone --depth 1 --branch v0.18.1 https://github.com/spack/spack.git
    cd spack
    echo "export SPACK_ROOT=$PWD" >> $HOME/.bashrc
    echo "export PATH=\$SPACK_ROOT/bin:\$PATH" >> $HOME/.bashrc
    # clone seissol-spack-aid and add the repository
    git clone --branch supermuc_NG https://github.com/SeisSol/seissol-spack-aid.git
    cd seissol-spack-aid
    spack repo add ./spack
    spack compiler find


Then update ``~/.spack/packages.yaml`` as follow:

.. code-block:: yaml

    packages:
      python:
        externals:
        - spec: python@3.10.1
          buildable: False
          modules:
           - python/3.10.1-cdl
      cmake:
        buildable: false
        externals:
        - spec: cmake@3.22.1
          modules:
           - cmake/3.22.1
      mpich:
        buildable: false
        externals:
        - spec: mpich@7.7.18
          modules:
          - cray-mpich/7.7.18
        - spec: mpich@7.7.20
          modules:
          - cray-mpich/7.7.20
      all:
        providers:
          mpi: [mpich]


Finally, install seissol-env with 

.. code-block:: bash

    spack install -j 8 seissol-env %gcc@11.2.0 ^mpich

and create a module with:

.. code-block:: bash

    spack module tcl refresh seissol-env@develop%%gcc@11.2.0

To access the module at start up, add to your ``~/.bashrc``:

.. code-block:: bash

    module use $SPACK_ROOT/share/spack/modules/cray-cnl7-ivybridge/

Finally, install SeisSol with cmake, as usual, with ``-DHOST_ARCH=hsw`` and ``-DCOMMTHREAD=OFF``.
