..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _boundary_conditions:

Boundary conditions
===================

Every face of the mesh carries a face type. Interior faces and dynamic rupture faces
join two cells; the remaining types close the domain and are evaluated from the cell
behind the face alone. Which mesh tag maps to which face type is described in
:doc:`puml-mesh-format`; this page describes what the types do.

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Face type
     - Default tag
     - What it does
   * - ``regular``
     - 0, 6
     - An interior face between two cells.
   * - ``freeSurface``
     - 1
     - The traction vanishes on the face.
   * - ``freeSurfaceGravity``
     - 2
     - The surface may move, with the hydrostatic pressure of its own elevation
       acting on it. Requires a cell without shear, see below.
   * - ``dynamicRupture``
     - 3, ≥ 65
     - A frictional interface; see :doc:`dynamic-rupture`.
   * - ``dirichlet``
     - 4
     - The exterior state is an affine function of the interior one, read from a
       file.
   * - ``outflow``
     - 5
     - Absorbing: outgoing waves leave, nothing comes back in.
   * - ``analytical``
     - 7
     - The exterior state is the analytical solution of the scenario.

Not every boundary condition is available in every build. A material model states
which face types are defined for it, and a solver states which ones its kernels
implement; a mesh that uses one that is not available is rejected at startup, with
the reason. The free surface with gravity additionally has to be checked per cell,
since its precondition is a property of the material in the cell rather than of the
material model.

Free surface
------------

The traction on the face is zero. Nothing has to be configured.

Free surface with gravity
-------------------------

The surface is free to move, and the hydrostatic pressure of its own elevation acts
on it. This is the ocean surface in a tsunami setup. The surface elevation follows an
ODE that is integrated along with the timestep, and its potential energy is reported
in the :doc:`energy-output`.

The ODE is closed with a single pressure and a scalar acoustic impedance, so the cell
behind the face must not carry shear. In an elastic run this means the cell needs a
vanishing shear modulus — which is how an ocean column is modelled, coupled to the
solid through the same solver. A face of this type on a cell with shear is rejected
at startup.

The gravitational acceleration is set in the equations block:

.. code-block:: Fortran

  &equations
  GravitationalAcceleration = 9.81
  /

Absorbing
---------

Only the outgoing characteristics leave the domain and nothing enters. This is exact
for a wave arriving perpendicular to the face and increasingly approximate the more
oblique the incidence, so place absorbing boundaries far enough from the region of
interest.

Dirichlet
---------

The state in the ghost cell behind the face is an affine function of the state inside:

.. math::

   q_\text{ghost} = A \, q_\text{inside} + b

Both :math:`A` and :math:`b` are read from an :doc:`easi` file, named in the equations
block:

.. code-block:: Fortran

  &equations
  BoundaryFileName = 'boundary.yaml'
  /

The condition is sampled once per face, at the face barycenter, so it may vary from
face to face but not within one.

Terms
~~~~~

The entries of :math:`A` are named ``map_{to}_{from}`` and those of :math:`b` are
named ``const_{to}``, where both names are quantities of the material at hand. For an
elastic material these are

.. code-block:: text

  s_xx  s_yy  s_zz  s_xy  s_yz  s_xz  v1  v2  v3

Entries that are not supplied default to the identity for :math:`A` and to zero for
:math:`b`, so a file only states what it changes. A term whose name is not a quantity
of the material is an error, and the message lists the names it could have been.

The offset :math:`b` prescribes a state over the timestep and is weighted with the
timestep width accordingly.

Frame
~~~~~

By default the condition is stated in global coordinates. A condition that mirrors or
fixes a direction is diagonal in the face-aligned basis instead, and on a boundary
that is not axis-aligned it cannot be written down globally at all. Setting

.. code-block:: yaml

  frame: 1

states the condition in the face-aligned basis, where the first axis is the face
normal, the second the strike direction and the third the dip direction. The
quantities keep their positions but refer to those axes, so for an elastic material

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Name
     - Means, with ``frame: 1``
   * - ``s_xx``
     - normal stress :math:`\sigma_{nn}`
   * - ``s_yy``, ``s_zz``
     - the in-plane normal stresses
   * - ``s_xy``, ``s_xz``
     - the shear tractions :math:`\sigma_{ns}`, :math:`\sigma_{nd}`
   * - ``s_yz``
     - the in-plane shear stress
   * - ``v1``
     - normal velocity
   * - ``v2``, ``v3``
     - the tangential velocities

``frame: 0``, the default, states the condition in global coordinates.

Examples
~~~~~~~~

A rigid wall, which reflects the normal velocity and lets the tangential motion slide:

.. code-block:: yaml

  !ConstantMap
    map:
      frame: 1
      map_v1_v1: -1.0

A symmetry plane, which additionally releases the shear tractions, so that half of a
symmetric domain can be meshed:

.. code-block:: yaml

  !ConstantMap
    map:
      frame: 1
      map_v1_v1: -1.0
      map_s_xy_s_xy: -1.0
      map_s_xz_s_xz: -1.0

A prescribed normal traction, stated globally on a boundary whose normal is the x
axis:

.. code-block:: yaml

  !ConstantMap
    map:
      const_s_xx: 1.0e6

Analytical
----------

The state behind the face is the analytical solution of the configured scenario,
evaluated at the face nodes and integrated over the timestep. It is used to measure
convergence against a known solution rather than to model a physical boundary.

The scenario has to be one that can be evaluated at an arbitrary time, not only at
:math:`t = 0`. Of the scenarios in :doc:`initial-condition` every analytically given
one qualifies; a scenario read from a file does not, since it supplies a state at
:math:`t = 0` alone. A mesh with analytical faces combined with such a scenario is
rejected at startup.
