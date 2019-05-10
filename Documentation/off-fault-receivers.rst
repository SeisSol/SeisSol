Off fault receivers
===================

Introduction
------------

Ascii receivers are enabled using the namelist output. Here is a
commented example:

.. code-block:: Fortran

  &Output
  pickdt = 0.01 ! Pickpoint Sampling
  pickDtType = 1 ! Pickpoint Type
  nRecordPoints = 22 ! number of Record points which are read from file
  RFileName = 'receivers.dat' ! Record Points in extra file
  iOutputMask = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
  /

If pickDtType = 2, output is generated every N time steps, where N is
set by pickdt. If pickDtType = 1, output is generated every pickdt
second.

receivers.dat is an ascii file describing the receivers coordinates in
the form:

::

  x1 y1 z1
  x2 y2 z2
  (...)
  xn yn zn

:ref:`wavefield-iouputmask` acts as well on the wavefield output.

Placing free-surface receivers
------------------------------

Placing receivers on the free-surface requires special care when a
realistic topography is used. The procedure to move receivers exactly to
the surface is described
`here <https://github.com/SeisSol/Meshing/tree/master/place_receivers>`__.

compiling place_receivers on supermuc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   module load netcdf/mpi
   export PKG_CONFIG_PATH=$NETCDF_BASE/lib/pkgconfig/:$PKG_CONFIG_PATH
   scons 
