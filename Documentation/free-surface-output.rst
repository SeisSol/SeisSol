Free surface output
===================

Introduction
------------

Velocities and ground deformations can be imaged as a surface
representation, in a file that can be opened in paraview. Threads or
nodes can be dedicated to write this output (see :ref:`asynchronous-output`),
but it is usually not necessary. The output to enabled in the Output
namelist:

.. code-block:: Fortran

  &Output
  SurfaceOutput = 1
  SurfaceOutputRefinement = 1
  SurfaceOutputInterval = 0.5
  /

If ``SurfaceOutputRefinement = 0``, one triangle is outputted for each
mesh cell. The unknowns are evaluated at the center of each cell.
``SurfaceOutputRefinement = 1`` subdivides each triangle, into 4
subtriangles. Higher SurfaceOutputRefinement would further subdivide
each subtriangle.

variables
---------

   | **u**, **v**, **w**: ground velocities, x y and z components
   | **U**, **V**, **W**: ground displacements, x y and z components
