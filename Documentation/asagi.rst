ASAGI
=====

The software package `ASAGI <https://github.com/TUM-I5/ASAGI>`__ can be
used to map gridded material properties of the domain to the mesh used
for a SeisSol simulation. As input a netcdf file with x-y-z coordinates
and corresponding material properties :math:`(\rho,\mu, \lambda )` is
required. By using
`asagiconv <https://github.com/SeisSol/SeisSol/tree/master/preprocessing/science/asagiconv>`__
in a pre-processing step it is also possible to create a netcdf file
containing the velocity model directly from the Community Velocity Model
(CVM-H).


.. _installing_ASAGI:

Installing ASAGI
----------------

Be careful that the python and gcc package is the same as for the
compilation of SeisSol in a later step!

example on SuperMuc
~~~~~~~~~~~~~~~~~~~

-  load the following modules (order matters!)

.. code-block:: bash

   module load mpi.ibm 
   module load netcdf/mpi
   module load hdf5/mpi/1.8.18
   module load intel/17.0
   module load gcc
   module load cmake

-  get the repository

On a cluster without restricted access to outside sources, you could then clone the repository using:

.. code-block:: bash

   git clone --recursive https://github.com/TUM-I5/ASAGI.github
   
On supermuc, you have to setup port forwarding as described in :ref:`compile_run_supermuc`.

Then you can clone the project with 

.. code-block:: bash

   git clone git@github.com:TUM-I5/ASAGI.git

followed by running ``fix_submodules`` to clone and setup the submodules.

-  set compiler options:

.. code-block:: bash

   export FC=mpif90
   export CXX=mpiCC
   export CC=mpicc

-  install:

.. code-block:: bash

   mkdir build
   cd build
   export CMAKE_PREFIX_PATH=$NETCDF_BASE
   cmake ../ -DCMAKE_INSTALL_PREFIX=<path_to_ASAGI>/build
   make -j8
   make install

-  set the following paths

.. code-block:: bash

   export PKG_CONFIG_PATH=<path_to_ASAGI>/build/lib/pkgconfig
   export LD_LIBRARY_PATH=<path_to_ASAGI>/build/lib

building SeisSol with ASAGI support
-----------------------------------

Simply add the following lines to the scons parameter file and make sure
you use the same python and gcc package as for the compilation with
ASAGI.

.. code-block:: bash

   asagi=yes
   zlibDir=<path_to_ASAGI>/build/lib/

**Known issues:** “can not find Asagi” while compiling SeisSol

There are a couple of options that can be checked:

-  Is SeisSol compiled with a different python package?
-  Are the paths to ASAGI correctly included? Check
   ``echo $PKG_CONFIG_PATH`` and ``echo $LD_LIBRARY_PATH``
-  When re-installing ASAGI again it might also help to remove the
   temporary files .sconf_temp/ and .sconsign.dblite within the SeisSol
   folder

generating the netcdf input file
--------------------------------

using asagiconv
~~~~~~~~~~~~~~~

Asagiconv (Located
`here <https://github.com/SeisSol/SeisSol/tree/master/preprocessing/science/asagiconv>`__)
allow querying data, vizualising and exporting to netcdf data from the
3D Velocity Model for Southern California. For more detail, see `ASAGI
docu <http://www.seissol.org/sites/default/files/asagi.pdf>`__.

velocity models given as structured grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Asagi expects a 3d structured grid netcdf file. Such a file can be
  generated from an ASCII file using the command:
  ``ncgen -b asagi_example.txt``
| Here is a typical example for the ASCII file:

::

   netcdf asagi_example {
   types:
     compound material {
       float rho ;
       float mu ;
       float lambda ;
     }; // material
   dimensions:
       x = 3 ; // Number of points in x-direction
       y = 2 ; // Number of points in y-direction
       z = 1 ; // Number of points in z-direction
   variables:
       float x(x) ;
       float y(y) ;
       float z(z);
       material data(z, y, x) ;
   data:
     x = 2, 2.5, 3 ; // Grid points in x-direction (must have the same spacing)
     y = -1, 0 ; // Grid points in y-direction (must have the same spacing)
     z = 0 ; // Grid points in z-direction (must have the same spacing)

     data =
     {1, -1, 10}, // rho,mu,lambda for x0, y0, z0
     {2, -2, 11}, // rho,mu,lambda for x1, y0, z0
     {3, -3, 12}, // rho,mu,lambda for x2, y0, z0
     {4, -4, 13}, // rho,mu,lambda for x0, y1, z0
     {5, -5, 14}, // rho,mu,lambda for x1, y1, z0
     {6, -6, 15} ; // rho,mu,lambda for x2, y1, z0
   }

Additionally, the netcdf file can be directly created using matlab or
python.

SeisSol parameter file
----------------------


A simple example file setting the elastic properties using EASI can be
found
`here <https://github.com/SeisSol/easi/blob/master/examples/101_asagi.yaml>`__.

Such a file would be called adding in the namelist equation:

.. code-block:: fortran

   MaterialFileName = 101_asagi.yaml

Further information
-------------------

For further information, the use of asagiconv and asagi and its
compilation, please see: `ASAGI
docu <http://www.seissol.org/sites/default/files/asagi.pdf>`__.

Known issues
------------

There is a bug, when using ASAGI with MPI. A workaround is described in https://github.com/SeisSol/SeisSol/issues/46.
