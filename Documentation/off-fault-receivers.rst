.. _off_fault_receivers:

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
  RFileName = 'receivers.dat' ! Record Points in extra file
  /

If pickDtType = 2, the output is generated every N time steps, where N is
set by pickdt. If pickDtType = 1, output is generated every pickdt
second.

receivers.dat is an ASCII file describing the coordinates of the receivers in
the form:

::

  x1 y1 z1
  x2 y2 z2
  (...)
  xn yn zn


The receivers files contain the time-histories of the stress tensor (6 variables) and the particle velocities (3).
Currently, there is no way to write only a subset of these variables.

Receiver time-histories are normally written at the end of the simulation, but this can be changed by setting :code:`ReceiverOutputInterval = 10.0`. Here the output is written every ReceiverOutputInterval seconds.

Rotational Output
-----------------
You can additionally choose to write the rotation of the velocity field by setting :code:`ReceiverComputeRotation=1` in the parameter file.
The rotation of the vector field is defined as :math:`\text{rot} v = \begin{pmatrix} \partial_2 v_3 - \partial_3 v_2 \\ \partial_3 v_1 - \partial_1 v_3 \\ \partial_1 v_2 - \partial_2 v_1 \\ \end{pmatrix}`.


Placing free-surface receivers
------------------------------

Placing receivers on the free-surface requires special care when a
realistic topography is used. The procedure to move receivers exactly to
the surface is described
`here <https://github.com/SeisSol/Meshing/tree/master/place_receivers>`__.

Compiling place_receivers on SuperMUC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Load the relevant :ref:`modules <compile_run_supermuc>`.

.. code-block:: bash
  
  git submodule update --init
  mkdir build && cd build
  cmake ..
  make -j

