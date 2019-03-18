Wave field output
=================

Introduction
------------

The wavefield can be written in a hdf5 file, in order to visualize it in
paraview. To speed-up the process, it is recommanded to dedicate a few
nodes to the writing tasks (see :ref:`asynchronous-output`).

Refinement
----------

| 0 (default): Refinement is disabled, i.e. only one cell is outputted
  for each element.
| 1: Refinement strategy is Face Extraction : 4 subcells per cell
| 2: Refinement strategy is Equal Face Area : 8 subcells per cell
| 3: Refinement strategy is Equal Face Area and Face Extraction : 32
  subcells per cell
| The unknowns are always evaluated at the centre of the subcell.

.. _wavefield-iouputmask:

iOutputMask
-----------

iOutputMask allow visualizing only part of the unknown. The stress
tensor (6), the velocities (3), the plastic strain tensor (6) and the
accumulated plastic strain eta, can be switched off or on by changing
the corresponding bit in the iOutputMask array.

OutputRegionBounds
------------------

Using the OutputRegionBounds parameter, under the &Output heading, in
the parameter.par file the user can define the region for which the
output is to be written. This region is provided in the following
format:

.. code-block:: Fortran

   OutputRegionBounds = xMin xMax yMin yMax zMin zMax

Example
-------

| Here is an example of wavefield output parametrisation:

.. code-block:: Fortran

   &Output
   OutputFile = '/output/prefix'
   iOutputMask = 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1
   OutputRegionBounds = -5e3 5e3 -10e3 10e3 -8e3 0e0
   Format = 6                          ! Format (6=hdf5, 10= no output)
   TimeInterval = 5.0                  ! Index of printed info at time
   printIntervalCriterion = 2          ! Criterion for index of printed info: 1=timesteps,2=time,3=timesteps+time
   refinement = 1
   /
