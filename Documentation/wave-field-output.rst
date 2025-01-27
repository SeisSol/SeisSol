..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

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
| The unknowns are always evaluated at the center of the subcell.

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

.. only:: comment

    IntegrationMask
    ---------------

    **Warning**: broken -> currently not compiling.
    To use this output, enable option INTEGRATE_QUANTITIES in cmake.
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
   /

High-Order VTKHDF Output
------------------------

The high-order wavefield output can be enabled by setting ``wavefieldvtkorder`` in the ``output`` section to a positive value, corresponding to the order of the output polynomial per cell.

