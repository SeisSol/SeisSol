Fault output
============

Introduction
------------

Two types of outputs are available for imaging the rupture. The rupture
characteristics can be assessed at a set of locations, using ASCII
receiver files, or overall the whole fault using files that can be
opened in ParaView. Threads or nodes can be dedicated to write this
latter output (see :ref:`asynchronous-output`),
but it is usually not necessary. The type of output generated depends on the value
of the variable `OutputPointType` of the DynamicRupture namelist:

.. code-block:: Fortran

  &DynamicRupture
  OutputPointType = 4
  /

| 0 : no output
| 3 : ASCII fault receivers
| 4 : paraview file
| 5 : both

Paraview output
---------------

This output is parametrized by the Elementwise namelist, example:

.. code-block:: Fortran

  &Elementwise
  printIntervalCriterion = 2 ! 1=iteration, 2=time
  printtimeinterval_sec = 1.0
  OutputMask = 1 1 1 1 1 1 1 1 1 1 1 1 !described herafter
  refinement_strategy = 1 ! or 2
  refinement = 1
  /

printIntervalCriterion
~~~~~~~~~~~~~~~~~~~~~~

If printIntervalCriterion = 1, the output is generated every N time steps,
where N is set by printInterval. This option only works with Global time
stepping. If printIntervalCriterion = 2, output is generated every
printtimeinterval_sec.

refinement
~~~~~~~~~~

If refinement = 0, one triangle is outputted for each mesh cell. The
unknowns are evaluated at the center of each cell. refinement = 1
subdivides each triangle, into 3 or 4 subtriangles depending on the
refinement_strategy. if refinement_strategy=1 splits each triangle into
3 triangles, sharing the triangle barycenter as a node. if
refinement_strategy=2, triangles are split into 4 triangles. Higher
refinement would further subdivide each subtriangle.

OutputMask
~~~~~~~~~~~

OutputMask allows visualizing only part of the unknown. The unknown can
be switched off or on by changing the corresponding bit in the
OutputMask array.

1. **SRs** and **SRd**: slip rates in strike and dip direction
2. **T_s**, **T_d**: transient shear stress in strike and dip
   direction, **P_n**: transient normal stress
3. **u_n**: normal velocity (note that there is no fault opening in SeisSol)
4. **Mud**: current friction, **StV**: state variable in case of RS friction
5. **Ts0**,\ **Td0**,\ **Pn0**: total stress, including initial stress
6. **Sls** and **Sld**: slip in strike and dip direction
7. **Vr**: rupture velocity, computed from the spatial derivatives
   of the rupture time
8. **ASl**: absolute slip
9. **PSR**: peak slip rate
10. **RT**: rupture time
11. **DS**: only with LSW, time at which ASl>D_c
12. **P_f** and **Tmp**: pore pressure and temperature

Ascii fault receivers
---------------------

The output is parametrized by the Pickpoint namelist, example:

.. code-block:: Fortran

  &Pickpoint
  printtimeinterval = 1
  OutputMask = 1 1 1 1 1 1 1 1 1 1 1 1 !described herafter
  nOutpoints = 24
  PPFileName = 'fault_receivers.dat'
  /

printtimeinterval
~~~~~~~~~~~~~~~~~

The output is generated every printtimeinterval (local) time step. Using
this output with local time-stepping may result in differently sampled
receiver files.

.. _ioutputmask-1:

iOutputMask
~~~~~~~~~~~

same as for ParaView output.


seissolxdmf
~~~~~~~~~~~

SeisSol paraview files (XDMF/Hdf5 or XDMF/binary files, describing the fault outputs and the free-surface/volume wavefield) can also be read using our python module `seissolxdmf <https://pypi.org/project/seissolxdmf/>`__.

Additional Ascii output
-----------------------

The rupture front can be outputted at every gauss points by enabling RF_output_on.

.. code-block:: Fortran

  &DynamicRupture
  RF_output_on = 0
  /

We nevertheless recommand using the Paraview fault output for vizualizing the rupture time instead of this ASCII output.
