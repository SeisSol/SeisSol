..
  SPDX-FileCopyrightText: 2018 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _wave_field_output:

Wave field output
=================

Introduction
------------

The wavefield can be written in an hdf5 file, in order to visualize it in
ParaView. To speed up the process, it is recommended to dedicate a few
nodes to the writing tasks (see :ref:`asynchronous-output`).

Refinement
----------

| 0 (default): Refinement is disabled, i.e. only one cell is outputted
  for each element.
| 1: Refinement strategy is Face Extraction: 4 subcells per cell
| 2: Refinement strategy is Equal Face Area: 8 subcells per cell
| 3: Refinement strategy is Equal Face Area and Face Extraction: 32
  subcells per cell
| By default, the unknowns are evaluated at the center of the subcell; see
  ``wavefieldprojection`` below for the alternative.

.. note::

   Up to and including SeisSol v1.3, the subcells of a refined wavefield output were sampled in a
   different vertex labeling than the one the output mesh was built with. As a result, the value
   written for a subcell was the solution at the center of one of its siblings -- for
   ``refinement = 1``, three of the four subcells of every element were affected, and for
   ``refinement = 2`` and ``3`` the inner subcells were sampled at a point that is not the center
   of any subcell at all. Only ``refinement = 0`` was unaffected, since the center of the whole
   element is invariant under that relabeling. Output written with ``refinement > 0`` by an older
   version therefore differs from what SeisSol produces now, and the difference is not a
   regression.

wavefieldprojection
-------------------

Controls how the solution is transferred onto the output points:

| ``pointwise`` (default): the solution is evaluated at the output points. For
  ``wavefieldvtkorder = -1`` these are the subcell centers, so this reproduces the classic
  wavefield output.
| ``l2``: the solution is projected onto the output space in the L2 sense. For
  ``wavefieldvtkorder = -1`` this is the average over each subcell, which is conservative but
  differs from the point value by O(h^2).


.. _wavefield-iouputmask:

iOutputMask
-----------

iOutputMask allows switching on and off the writing of SeisSol unknowns.
The 6 first digits controls the components of the stress tensor
(sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, and sigma_xz),
and the 3 last digits the velocity components (u, v, w).
When using poroelasticity, 4 more flags are added, for pore pressure (p) and fluid velocities (u_f, v_f, w_f).

iPlasticityMask
---------------

iPlasticityMask allows switching on and off the writing of plasticity variables.
The 6 first digits controls the components of the off-fault plastic
strain tensor (ep_xx, ep_yy, ep_zz, ep_xy, ep_yz, and ep_xz),
and the last one the accumulated plastic strain (eta).

IntegrationMask
---------------

IntegrationMask allows switching on and off the writing of time integrated SeisSol unknowns.
The 6 first digits control the components of the time integrated stress tensor
(int_sigma_xx, int_sigma_yy, int_sigma_zz, int_sigma_xy, int_sigma_yz, and int_sigma_xz),
and the 3 last digits the displacement components (displacement_x, displacement_y, displacement_z).
Note that this output is associated with the prefix-low.xdmf file, and can only output
the cell average quantities.


OutputRegionBounds
------------------

Using the OutputRegionBounds parameter, under the &Output heading, in
the parameter.par file, the user can define the region for which the
output is to be written. This region is provided in the following
format:

.. code-block:: Fortran

   OutputRegionBounds = xMin xMax yMin yMax zMin zMax

OutputGroups
------------------

Similar to the previous parameter, OutputGroups can be used to whitelist a set of
mesh groups (as specified in the xdmf mesh file) that are included in the wavefield output.
Cells whose group is not mentioned are not included in the output.
This feature works with OutputRegionBounds, only cells that satisfy both criteria are included.
It looks like this:

.. code-block:: Fortran

   OutputGroups = 1 2 ! only include groups 1 and 2

Example
-------

| Here is an example of wavefield output parametrization:

.. code-block:: Fortran

   &Output
   OutputFile = '/output/prefix'
   iOutputMask     = 0 0 0 0 0 0 1 1 1
   iPlasticityMask = 0 0 0 0 0 0 1
   OutputRegionBounds = -5e3 5e3 -10e3 10e3 -8e3 0e0
   Format = 6                          ! Format (6=hdf5, 10= no output)
   TimeInterval = 5.0                  ! Index of printed info at time
   printIntervalCriterion = 2          ! Criterion for index of printed info: 1=timesteps,2=time,3=timesteps+time
   refinement = 1
   wavefieldvtkorder = -1
   wavefieldprojection = 'pointwise'
   wavefieldtimeseries = 'snapshot'
   /

File groupings
--------------

``wavefieldtimeseries`` overrides ``outputtimeseries`` for the wavefield output,
so that it can be written as one file for the whole run while the other outputs
stay one file per step, or the other way round. It takes effect only with
``wavefieldvtkorder`` set; see :ref:`io_infrastructure` for what the groupings
are.

At ``wavefieldvtkorder = -1`` and ``0`` the corners a cell shares with its
neighbors are written once, which makes the point array as large as the mesh
rather than as large as the mesh times the number of cells a vertex touches. The
merging happens within a rank, and can be switched off with
``SEISSOL_IO_VERTEXFILTER=0``.

High-Order VTKHDF Output
------------------------

The high-order wavefield output can be enabled by setting ``wavefieldvtkorder`` in the ``output`` section to a positive value, corresponding to the order of the output polynomial per cell.

Derived outputs
---------------

``wavefieldscript`` names a program whose outputs are written along with the wavefield, at the
same points. It is an sderiv module (``sderiv:file`` or a ``.sderiv`` file) or a Lua module
(``lua:file`` or a ``.lua`` file), written pointwise. A program reads the quantities by name:

| ``q`` -- a quantity of the solution, e.g. ``v1`` or ``s_xx``
| ``q_r0``, ``q_r1``, ``q_r2`` -- its derivative along the reference coordinates of the cell
| ``dx_q``, ``dy_q``, ``dz_q`` -- its derivative in space
| ``int_q`` -- the time integral of ``q`` (and its derivatives as above, e.g. ``dx_int_v1``)
| ``ep_xx`` ... ``eta`` -- the plastic strain, with plasticity
| ``jinv00`` ... ``jinv22`` -- the inverse Jacobian of the cell, ``jinvkd`` = d xi_k / d x_d
| ``x``, ``y``, ``z``, ``t`` -- the output point and the time
| ``dt`` -- the time since the previous evaluation of the point

The values and derivatives are taken from the coefficients of the cell directly, so a program
pays for each of them once however many outputs read it. The built-in outputs (the quantities,
``int-`` quantities, strain, rotation and plastic strain) are computed the same way.

A program without state is evaluated when the output is written. A program with state follows
every time step of the cells (on the CPU; on a GPU, it is evaluated when written), so that e.g. a
maximum over time does not miss what happens between two outputs:

.. code-block:: text

   # pgv.sderiv
   state pgv = 0.0
   out def pgv = max(pgv, sqrt(v1*v1 + v2*v2 + v3*v3))
   # the displacement, integrated over the time steps
   state u1 = 0.0
   out def u1 = u1 + v1 * dt
   out def divv = dx_v1 + dy_v2 + dz_v3

The same in Lua; a returned table names the outputs, and ``M.state`` declares the state:

.. code-block:: lua

   local M = {}
   M.state = { pgv = 0.0 }
   function M.evaluate(fields, v1, v2, v3, pgv)
     return { pgv = math.max(pgv, math.sqrt(v1*v1 + v2*v2 + v3*v3)) }
   end
   return M

.. code-block:: Fortran

   &Output
   wavefieldscript = 'sderiv:pgv.sderiv'
   /

A state is not written to checkpoints; a restarted run starts it over.
