..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _free_surface_output:

Free surface output
===================

Introduction
------------

Velocities and ground deformations can be imaged as a surface
representation, in a file that can be opened in ParaView. Threads or
nodes can be dedicated to write this output (see :ref:`asynchronous-output`),
but it is usually not necessary. The output to enabled in the Output
namelist:

.. code-block:: Fortran

  &Output
  SurfaceOutput = 1
  SurfaceOutputRefinement = 1
  SurfaceOutputInterval = 0.5
  surfacevtkorder = -1
  /

If ``SurfaceOutputRefinement = 0``, one triangle is outputted for each
mesh cell. The unknowns are evaluated at the center of each cell.
``SurfaceOutputRefinement = 1`` subdivides each triangle, into 4
subtriangles. Higher SurfaceOutputRefinement would further subdivide
each subtriangle.

variables
---------

   | **v1**, **v2**, **v3**: ground velocities, x y and z components
   | **u1**, **u2**, **u3**: ground displacements, x y and z components

Additionally, the writer outputs a quantity called "locationFlag", which has the values
0 and 1 when at the elastic or acoustic side of an elastic-acoustic interface.
In this way, we can distinguish between both sides of the interface even though they have the same coordinates.
It has the value 2 for an ordinary free surface boundary condition and the value 3 for a free surface with gravity
boundary condition.
This value can be used to filter the output (which contains all these surfaces), for example using Paraview's Threshold filter.

High-Order VTKHDF Output
------------------------

The high-order free surface output can be enabled by setting ``surfacevtkorder`` in the ``output`` section to a positive value, corresponding to the order of the output polynomial per cell.

