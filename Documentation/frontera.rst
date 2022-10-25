.. _compile_run_frontera:


Frontera
========


Add these lines to ``~/.bashrc``:

::

    module switch python3 python3/3.9.2
    module use /work2/09160/ulrich/frontera/spack/share/spack/modules/linux-centos7-cascadelake
    module load seissol-env-develop-intel-19.1.1.217-x52n3zf
    export CC=mpiicc
    export CXX=mpiicpc
    export FC=mpiifort

This will load a preinstalled seissol-env module.

Alternatively (and for reference), to compile seissol-env on Frontera, follow the procedure below:

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

Following the workaround proposed in https://github.com/spack/spack/issues/10308, precise the module of the intel compilers in ``~/.spack/linux/compilers.yaml`` by changing ``modules: []`` to ``modules: ['intel/19.1.1']``.

Then, update ``~/.spack/packages.yaml`` as follow:

.. code-block:: yaml

    packages:
      autoconf:
        externals:
        - spec: autoconf@2.69
          prefix: /opt/apps/autotools/1.2
      python:
        externals:
        - spec: python@3.9.2
          buildable: false
          prefix: /opt/apps/intel19/python3/3.9.2
          modules:
          - python3/3.9.2
      intel-mpi:
        buildable: false
        externals:
        - spec: intel-mpi@2019.0.9
          modules:
          - impi/19.0.9
      all:
        providers:
          mpi: [intel-mpi]

(note that the compilation was not successful with trying to add the cmake/3.24.2 module to packages.yaml).

Finally, install seissol-env with 

.. code-block:: bash

    spack install -j 16 seissol-env %intel@19.1.1.217 ^intel-mpi

and create a module with:

.. code-block:: bash

    spack module tcl refresh seissol-env

To access the module at start up, add to your ``~/.bashrc``:

.. code-block:: bash

    module use $SPACK_ROOT/share/spack/modules/linux-centos7-cascadelake/

Finally, install SeisSol with cmake, as usual, with ``-DHOST_ARCH=skx``.
