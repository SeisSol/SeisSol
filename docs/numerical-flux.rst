..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _numerical_flux:

Numerical flux
==============

At a face between two cells, the discontinuous Galerkin method sees two values of the
solution, one from each cell. The numerical flux turns them into one flux through the face.
SeisSol takes it from a state on the face that both cells share, the Riemann state. Which
flux to use is set in the ``&Equations`` section of the parameter file, for the faces next to
the fault separately:

.. code-block:: Fortran

  &Equations
  numflux = 'godunov'
  numfluxnearfault = 'godunov'
  /

``numflux`` applies to all regular faces. ``numfluxnearfault`` applies to the regular faces
of the cells that have a dynamic rupture face: the faces between two such cells and those
between such a cell and its other neighbors. The dynamic rupture faces themselves take their
state from the friction law.

Fluxes
------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Value
     - Flux
   * - ``'godunov'``, ``'upwind'``
     - The exact solution of the Riemann problem between both cells (the default): the
       impedance-weighted means of their traction and velocity, and an upwind damping of the
       jumps between them.
   * - ``'centered'``
     - The mean of both sides, without damping.
   * - ``'rusanov'``
     - The mean, and a damping of the jumps of all quantities with the largest wave speed of
       both sides (the local Lax-Friedrichs flux).
   * - ``'mixed'``
     - ``fluxgpart`` times the Godunov flux and the rest the centered one, and ``fluxrpart``
       times the damping of the Rusanov flux; near the fault ``fluxgpartnearfault`` and
       ``fluxrpartnearfault``.
   * - ``'godunov-shearcentered'``
     - The Godunov flux without the damping of the tangential traction jumps
       (``fluxdamping = '1 1 1 0'``), and a tenth of the damping of the Rusanov flux; meant for
       ``numfluxnearfault``, see below. It sets these parts itself; to change them, write them
       out with ``'mixed'``.

Every value but ``'mixed'`` and ``'godunov-shearcentered'`` ignores ``fluxgpart`` and
``fluxrpart``. The log shows at startup which flux each kind of face takes, with its parts.

The shares of the damping
-------------------------

The damping of the Godunov flux has two parts: the traction of the Riemann state contains a
damping of the jumps of the velocity, and its velocity a damping of the jumps of the traction.
``fluxdamping`` (and ``fluxdampingnearfault`` for the faces near the fault) scales them, each
in its normal and its tangential direction:

.. code-block:: Fortran

  fluxdamping = 'vn vt tn tt'

``vn`` and ``vt`` are the shares of the jumps of the normal and the tangential velocity,
``tn`` and ``tt`` those of the normal and the tangential traction. A share of 1 keeps the
damping of the Godunov flux (the default is ``'1 1 1 1'``), 0 removes it, and values in
between or above are allowed. On a face within one material, ``'0 0 0 0'`` is the centered
flux; at a material contrast, the impedance-weighted means remain.

* A poroelastic material takes two more shares, ``'vn vt tn tt vf p'``, for the jumps of the
  normal velocity of the fluid and of its pressure. Given four shares, the fluid takes those of
  the normal jumps of the solid, ``vn`` and ``tn``.
* An acoustic material has only the normal shares ``vn`` and ``tn``.
* Where the material couples components of different directions (anisotropy) or the solid with
  the fluid, the damping between two components is scaled by the geometric mean of their
  shares. That keeps it a damping: with shares of at least 0, no face creates energy.

With ``'mixed'``, the shares apply to its Godunov part; ``'centered'`` and ``'rusanov'`` have
none.

.. _numerical_flux_face_map:

The flux by face tags and mesh groups
-------------------------------------

The face map (``facemapfile``, see :doc:`puml-mesh-format`) can give regular faces a flux of
their own. An entry of ``regular`` that is a map selects faces by their ``tags``, by the mesh
``groups`` of their cells, or by both, and gives their flux with the keys of the parameter file:

.. code-block:: yaml

   regular:
     - 0, 6
     - tags: 100-110             # regular faces of their own tags, with their own flux
       numflux: godunov-shearcentered
     - groups: 2, 5-7            # the regular faces of the cells of these mesh groups
       numfluxnearfault: mixed
       fluxgpartnearfault: 0.5
   freeSurface: 1
   dynamicRupture: 3, 65-
   outflow: 5

* A face is selected by ``groups`` if the cell on either side of it belongs to one of them; with
  ``tags`` and ``groups``, it has to match both.
* ``numflux`` and its parts set the flux of the selected faces away from the fault,
  ``numfluxnearfault`` and its parts that of those near it; an entry that names only one of them
  leaves the other faces as they are. The rules of the parameter file apply as well: a preset
  sets its parts itself, and only ``'mixed'`` takes ``fluxgpart`` and ``fluxrpart``.
* The first entry that selects a face and sets a flux for it decides; the faces that no entry
  selects take the flux of the parameter file.
* Both cells of a face see its tag and both mesh groups, and take the same flux.
* Only regular faces take a flux of their own; the other face types keep the Godunov flux (see
  below).

The log lists at startup how many regular faces take which flux, and from where.

Faces with the Godunov flux
---------------------------

Some faces keep the Godunov flux whatever the parameters say, since their Riemann state
imposes a condition there:

* all boundary faces (free surface, free surface with gravity, Dirichlet, outflow, analytical;
  see :doc:`boundary-conditions`), where the state is that of the boundary condition: a smaller
  damping would, for instance, let an outflow boundary reflect;
* faces between an elastic and an acoustic cell, where the state keeps the shear traction
  zero;
* faces between cells of different configurations, and between poroelastic cells with a
  hydraulic skin between them (see :doc:`physical-models`), whose states couple the two sides
  in their own way.

Stability
---------

With the Godunov flux and any shares of at least 0, and with the centered flux, no face creates
energy; neither does the damping of the Rusanov flux within one material. The time stepping,
however, needs the damping as well: a face that leaves some jumps undamped can make the scheme
unstable, depending on the order. On the TPV5 benchmark (local time stepping, CFL 0.5), with
these fluxes on the faces near the fault:

* ``'centered'`` grew without bound already with order 4, also with CFL 0.1.
* Without the damping of the tangential traction jumps (``fluxdampingnearfault = '1 1 1 0'``),
  orders 3 and 4 stayed stable, while orders 5 and 6 grew without bound within one to two
  seconds.
* A tenth of the damping of the Rusanov flux (``fluxrpartnearfault = 0.1`` with
  ``fluxgpartnearfault = 1``) kept orders 4 to 6 stable over 8 s, which is what
  ``'godunov-shearcentered'`` adds. Keeping a part of the damping of the tangential traction
  jumps instead (0.1 with order 5, 0.2 with order 6) also suffices, but brings back much of the
  drift described below.

Linear slip weakening on meshes that are not symmetric about the fault
----------------------------------------------------------------------

With linear slip weakening, the strength of the fault is its friction coefficient times the
normal stress. Where the cells on the two sides of the fault are not mirror images of each
other, as in most unstructured meshes, the normal stress on the fault picks up a spurious part
during the rupture, which changes the strength, which changes the rupture. On faces next to the
fault that are oblique to it, the damping of the tangential traction jumps is the main source
of this spurious part. It depends on the polynomial order and the mesh, and shows, for instance,
in a seismic moment that changes with them more than the resolution would explain.

``numfluxnearfault = 'godunov-shearcentered'`` reduces the drift. On TPV5 with smooth
patches, the seismic moment after 8 s deviated from that of a reference with the normal stress
ten times larger and the friction coefficients ten times smaller (the same strength, ten times
less sensitive to the normal stress) by:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Order
     - ``'godunov'``
     - ``'godunov-shearcentered'``
   * - 4
     - +13.7 %
     - +0.2 %
   * - 5
     - +11.1 %
     - +0.9 %
   * - 6
     - +7.4 %
     - +1.4 %

The higher the order, the less it helps. It does not remove the cause either: the drift
returns where the normal stress matters more (twice the sensitivity of the benchmark already
lets it grow again). A mesh whose first layer of cells is symmetric about the fault removes it.
