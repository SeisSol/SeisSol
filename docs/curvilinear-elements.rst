..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Curvilinear elements
====================

A build with ``CURVILINEAR=ON`` lets a cell be curved. This page says what the
scheme does with a curved cell, where that needs data per point -- per point
the operator is formed at, per node of a face -- and where it does not, which
parts of the code had to learn about it, and what is not done yet. The
measurements behind the choices are summarised here; the commit messages of the
series carry the details.

Building and feeding it
-----------------------

``CURVILINEAR=ON`` needs ``MATERIAL_NODAL=ON`` and
``MATERIAL_OPERATOR=assembled``, runs on the CPU, and takes the media whose flux
decomposes into scalars per face node: elastic, acoustic and the attenuating
variants of both, not the anisotropic or the poroelastic one. cmake and the code
generator refuse the rest. A material that does not vary is carried as samples
that agree.

Which cells are curved, the mesh decides. PUMGen writes a mesh of gmsh whose
cells are of a higher order (from its branch ``davschneller/more-modernize`` on)
as it writes any other mesh -- the vertices in ``/geometry``, the cells by their
vertices in ``/connect`` -- and adds the other nodes of every cell:
``/geometry_ho``, in the order gmsh numbers them (its attribute
``node-ordering`` says ``gmsh``), ``/geometry_ho_offsets``, where those of every
cell start, and ``/order``, the order of every cell. The PUML reader reads them
as a cell array of its own length per cell, which PUML moves with its cell
through the partitioning (from its master on, which the series moves to).

Before the partitioning, the time step of a cell comes from its curved map
(``CellToVertexArray::fromPUML``). After it, the reader puts the nodes into the
vertex order it keeps and onto the lattice of the highest order of the mesh --
a cell of a lower order is the same map written with more nodes -- and hands
them over through ``MeshReader::setCurvedGeometry``: the nodes of every local
cell on the equispaced lattice of the geometry order
(``IsoparametricTransform::latticeNodes``; order two adds the six edge midpoints
to the vertices), the vertices first. The call checks that the vertices are the
cell's and that the cell is not turned inside out anywhere on a lattice of twice
its order.

A build without ``CURVILINEAR`` takes such a file only where all cells are
straight-sided -- their nodes within :math:`10^{-8}` of their longest edge of
where the straight-sided cell has them -- and then leaves the nodes aside; a
mesh with a curved cell it refuses before anything is partitioned.

On a straight-sided mesh a ``CURVILINEAR`` build computes what the build without
it computes: a plane wave through a cube gives the same error norms to
:math:`10^{-12}` and the same energies to :math:`10^{-15}`.

At order 4 a cell carries 9.0 kB of metric at the 125 operator points and 14.4 kB
of face rotations at the 10 nodes of each face, against 72 B and 1.4 kB for a
straight-sided cell.

The scheme on a curved cell
---------------------------

A cell is a map :math:`x(\xi)` from the reference tetrahedron, its Jacobian
:math:`J = \partial x / \partial \xi` varies inside it, and so does the metric
:math:`J^{-1}`, whose rows are the gradients of the reference coordinates. For
an isoparametric map of order :math:`g`, :math:`\det J` has the degree
:math:`3(g-1)` and the cofactor matrix :math:`\det(J)\,J^{-1}` the degree
:math:`2(g-1)`.

**The volume term is the strong form.** The predictor and the volume term of the
corrector both differentiate the field in reference coordinates, evaluate that
at the points the operator is formed at, apply the operator there --
:math:`\sum_j (J^{-1})_{ej}(\xi_n) A_j(\xi_n)`, the metric at the point together
with the material at the point, folded once per kernel -- and project back.
That is what a nodal material build did already for the predictor; the corrector
used the weak derivative instead, and with an operator that varies inside a cell
that is not consistent: the weak derivative is the derivative less the lift of
the cell's own trace, the operator inside the cell multiplies that lift, and the
face flux takes the trace off again with the operator at the face. The two only
cancel where they are one operator. Measured on one cell (scalar advection,
order 4, relative error of :math:`\partial_t q`):

.. list-table::
  :header-rows: 1
  :widths: 28 18 18 18 18

  * - :math:`h`
    - 0.4
    - 0.2
    - 0.1
    - 0.05
  * - curved, weak form
    - 7.7e+1
    - 7.1e+0
    - 5.9e+0
    - 5.5e+0
  * - curved, strong form
    - 1.0e-1
    - 9.4e-4
    - 1.2e-4
    - 1.5e-5
  * - varying material, weak form
    - 7.4e-1
    - 6.5e-1
    - 6.1e-1
    - 5.9e-1
  * - varying material, strong form
    - 4.4e-4
    - 4.3e-5
    - 4.7e-6
    - 5.4e-7

So the strong form is not a choice for curved cells, and it is the fix for the
nodal material as well: through the whole code, a constant state in a smoothly
varying medium drifted by 2 to 18 percent within 0.1 s with the weak form, and
by about as much on a mesh twice as fine; with the strong form it stays
constant to :math:`10^{-14}`. The local flux of every face now takes off the
normal flux of the cell's own trace (``toCorrectorForm``), a fault face
included, whose local flux is that subtraction alone.

**The mass matrix stays the one of the reference cell.** The scheme tests with
the basis functions divided by :math:`\det J`. A point source then is
:math:`\varphi(\xi_s) / \det J(\xi_s)`, the initial state is the projection of
its values at the quadrature points, and the face terms take on
:math:`|n| / \det J` at every face node, :math:`|n|` the surface Jacobian -- for
a straight-sided cell that is the flux scale :math:`2|S| / |J|` it always had.
Nothing that has the inverse mass matrix folded in -- ``kDivMT``, ``rDivM``,
``fMrT``, ``project2nFaceTo3m``, ``M3inv``, ``V3mTo2nTWDivM`` -- changes, and no
cell carries a dense mass matrix. The price is exact conservation, which the
scheme for a material that varies inside a cell does not have either, and
stability is not proven the way it is for the Galerkin scheme with
:math:`M_J`. Measured instead: the spectrum of the semi-discrete operator of 1D
acoustics on curved cells, with a varying material and the Godunov flux, has no
eigenvalue with a positive real part beyond round-off at orders 4 and 6, and in
3D a constant state stays constant on curved meshes.

**A curved face has its rotation and its scale per node.** Its normal turns
along it, so the rotation into face coordinates is one per node (``TNodes``),
taken from the frame the face has at the node, and so is the scale. The flux
scalars are formed per node as before.

**The time step** of a curved cell is the one of its vertices times its
relative thickness (``relativeThickness``): the smallest singular value of the
Jacobian over the one of the straight-sided cell, where it is smallest. A cell
squeezed to 5 percent of its straight Jacobian determinant blew up with any CFL
down to 0.1 before this.

**A curved fault has its frame per point.** Its normal is the one of the face
at the point, and its first tangent the one of the fault -- along the first edge
of the plane through its vertices -- turned into the plane the face has there
(``faultFrameAt``); both sides of a fault come to the same frame that way, on
one rank to the bit and across ranks to round-off. The fault rotates the state
of either side into the frame of every point and its imposed state back out of
it (``TinvTPoints``, ``TPoints``, which replace ``TinvT`` and ``T`` in a
``CURVILINEAR`` build) and lifts it with the scale of the side at the point,
:math:`\mp|n|/\det J`. Its initial stress and the slip it is made to have are
those of the frame of their point, and its parameters are queried at the points
on the curved face. Its output gives what it evaluates at an output point in
the frame there, and turns what it reads off the nearest quadrature point -- the
initial stress, the slip -- out of the frame of that point. The resampling is
not touched: it maps between the points of the reference face and does not know
the geometry.

**A curved fault subtracts the trace of a side itself.** In the strong form the
local flux of a fault face is the subtraction of the normal flux of the cell's
own trace, which it takes at the nodes of the face, while the fault lifts its
imposed state at its points. On a curved face the scale :math:`|n|/\det J` is
not a polynomial, the two quadratures differ, and a constant state across a
locked fault drifted by :math:`2.4 \cdot 10^{-5}` within 0.2 s. So in a
``CURVILINEAR`` build the fault lifts the imposed state less the side's own
trace, both integrated in time with the same weights, and the local flux leaves
fault faces alone (``FaultSubtractsOwnTrace``); the constant state then stays
constant to :math:`10^{-14}`, as across a plane fault.

Where it has to be per point
----------------------------

.. list-table::
  :header-rows: 1
  :widths: 22 14 64

  * - Part
    - Per point?
    - Why
  * - Face flux
    - yes, per node
    - Rotated, the operator of a curved face is not one matrix of the face
      times one matrix of the quantities, so neither the rotation nor the scale
      can be folded into what a face stores. This is the one place where it has
      to be pointwise in any formulation.
  * - Volume operator, predictor and corrector
    - yes, per operator point
    - The metric multiplies the derivative at the point. One could integrate it
      into per-cell matrices instead (``davschneller/config`` does), which are
      dense and 40.6 kB per cell at order 4; with the nodal operator already
      formed at points, the metric rides along in the assembly.
  * - Mass matrix
    - no
    - The reference mass matrix stays, see above. The Galerkin mass matrix
      :math:`M_J` would be dense per cell and folded into six matrix families.
  * - Point sources
    - one point
    - :math:`1 / \det J` at the source.
  * - Energies, error norms
    - per quadrature point
    - :math:`\det J` weighs every point.
  * - Material, initial state, receivers, output
    - through the map
    - Points are placed through the cell's map; a receiver finds its reference
      point by Newton's method, and derived output takes the metric at each
      point.
  * - Plasticity
    - no
    - It acts on values at points and does not differentiate.
  * - Dynamic rupture
    - yes, per point
    - The frame of a curved fault turns along it, so the rotation into it and
      out of it, the scale of the lift and the subtraction of the own trace are
      per point, see above, and so are the frames of the initial stress, of an
      imposed slip and of the output, and the weights of the fault energies.
  * - Boundary conditions with a ghost state
    - would be per node
    - Analytical, Dirichlet and gravity state a ghost cell and apply one matrix
      per face. Not done; refused on a curved mesh. Free surface and outflow are
      part of the flux and work.

What is not done
----------------

- Curved cells from the readers other than the one of PUML, and from a mesh of
  several kinds of cells (which PUMGen writes in the layout of VTKHDF).
- The boundary conditions with a ghost state on curved faces, see above.
- The fault output on a curved face: its points are on the face and take the
  frame of the quadrature point they read the fault at, but the triangles the
  elementwise output writes at order zero are plane between their corners.
- Ghost cells on other ranks are sampled straight-sided: their metadata carries
  their vertices only, so the material a face reads from a neighbour on
  another rank comes from sample points placed by the straight-sided map.
- The device path. A curved face's rotation per node has not been generated or
  checked with tensorforge.
- ``MATERIAL_OPERATOR=factored``, which would fold the metric at every
  application.
- Memory: the metric could be carried at the material samples and
  interpolated, as the material is, instead of at the operator points; the face
  rotation could be built in the kernel from a normal per node.

Checks
------

Unit checks, in either build where they apply: a straight-sided cell evaluated
at every point carries what it carries as one cell; the isoparametric map passes
through its nodes, its Jacobian is its derivative, its cofactor matrix is free
of divergence and gives the face normal (Nanson); on a curved cell the face
scale and frame per node agree with differencing the map; the volume term
differentiates a linear field exactly on a curved cell -- with the metric of one
point for the whole cell, a mode that has to vanish comes out at a sizeable
fraction of the derivative -- and a constant state stays constant; the strong
form with its face terms is the weak form for one material per cell, a fault
face included. On a fault: the lift of every point is the flux of the normal
there, with the rotation and the scale of the point; the frame of a plane face
is the one of the fault at every point, and the one of a curved face is the
normal there with the tangent of the fault turned into it; the fault subtracts
the trace of a side integrated in time; and the output finds the point of a
curved face nearest to a point it is given.

End to end, with a hook that is not part of the series: a cube of
:math:`n^3` hexahedra cut into six tetrahedra each, its vertices moved by the
smooth map :math:`x + \varepsilon \sin(\pi x)\sin(\pi y)\sin(\pi z)\,(1, 0.7,
-0.5)` with :math:`\varepsilon = 0.06`, which leaves its boundary where it is,
and either straight edges between them or the edge midpoints where the map puts
them, curving the cells. A Gaussian pulse runs for 0.15 s with outflow on the
boundary; four receivers inside are compared with a straight-sided mesh of
:math:`24^3` hexahedra, as the relative L2 difference of their time series
(largest over the receivers):

.. list-table::
  :header-rows: 1
  :widths: 20 20 20 20 20

  * - :math:`n`
    - straight edges, constant material
    - curved, constant material
    - straight edges, varying material
    - curved, varying material
  * - 4
    - 1.9e-2
    - 2.2e-2
    - 2.0e-2
    - 2.3e-2
  * - 8
    - 1.5e-3
    - 1.4e-3
    - 1.6e-3
    - 1.6e-3
  * - 12
    - 3.1e-4
    - 3.1e-4
    - 3.4e-4
    - 3.6e-4

Both converge at the rate of the order, and the curved cells are as accurate as
the straight ones. The mesh, the parameters and the hook are with the
measurement scripts that accompany the series.

The same cubes written as PUMGen writes cells of a higher order and read by the
PUML reader give the same differences to the reference (2.2e-2 for
:math:`n = 4`, 1.5e-3 for :math:`n = 8`); with the time step the reader takes for
them -- shorter by the relative thickness of the cells -- the series differ from
the ones of the hook by less than a fifth of that, and with one fixed time step
for both they agree to 1e-15.

**A curved boundary.** The unit ball meshed by gmsh at order two and converted by
PUMGen, its surface a free surface, oscillating in its fundamental torsional mode
:math:`{}_0T_2` (:math:`\mu = \rho = 1`, the frequency the root
:math:`x \approx 2.5011` of :math:`j_2(x) = x\, j_3(x)`) for two periods; five
receivers against the mode, as the relative L2 difference of their time series
(largest over the receivers), with the cells of order one -- the same vertices,
plane faces between them -- and of order two:

.. list-table::
  :header-rows: 1
  :widths: 25 25 25 25

  * - mesh size
    - cells
    - order one
    - order two
  * - 0.5
    - 256
    - 1.9e-1
    - 3.1e-3
  * - 0.35
    - 503
    - 1.2e-1
    - 1.3e-3
  * - 0.25
    - 1435
    - 5.1e-2
    - 2.7e-4

With plane faces the ball is the polyhedron inside the sphere, which is too
small, and the difference falls with the square of the mesh size (1.9 from the
coarsest mesh to the finest); with the cells of order two it falls faster than
with its cube (3.5).

**A curved fault.** The cubes again, in the layout PUMGen writes, with the faces
of a coordinate plane through the middle of the unmoved cube a fault, which the
map bends. Locked (cohesion :math:`-10^{10}`), in a periodic cube of
:math:`4^3` bent further (:math:`\varepsilon = 0.1`), a constant state stays
constant to :math:`2.5 \cdot 10^{-14}`
within 0.2 s, as across the same fault left plane (:math:`3.8 \cdot
10^{-14}`); with the subtraction of the own trace at the nodes of the face it
drifted by :math:`2.4 \cdot 10^{-5}`. The Gaussian pulse of the cubes above,
started on the fault (the plane :math:`z = 1/2` bent) and run across it, gives
the same differences to the reference as without the fault:

.. list-table::
  :header-rows: 1
  :widths: 25 25 25 25

  * - :math:`n`
    - without the fault
    - across the locked fault
    - between the two
  * - 4
    - 2.2e-2
    - 2.2e-2
    - 1.0e-3
  * - 8
    - 1.5e-3
    - 1.5e-3
    - 1.4e-4
  * - 12
    - 3.1e-4
    - 3.1e-4
    - 1.7e-5

With the plane :math:`x = 1/2` bent into the fault, so that strike and dip are
well defined, three receivers on it give the tractions; against the ones of the
stress of the reference at the same points, on the normal of the bent plane
there, the relative L2 difference is 3.9e-2, 1.8e-3 and 2.7e-4 for the normal
traction (1.5e-2, 1.4e-3, 1.2e-4 along strike; 1.2e-2, 9.8e-4, 1.7e-4 along
dip) for :math:`n = 4, 8, 12`. (The receivers on a fault evaluate the cells at
the start of the step whose end they write as the time; the comparison takes
that into account. It is so for a plane fault as well.) Under a uniform stress
the fault starts out under, the total tractions they give differ from the ones
on the normal of the bent plane by 7.1e-3, 2.0e-3 and 5.0e-4 -- what the faces
of order two make of its normal -- while the wavefield stays zero. On four ranks
the series are the ones of one rank, to the bit where ParMETIS keeps the fault
inside the partitions, and to :math:`10^{-15}` where the fault is the boundary
between two of them.

With slip the fault is made to have -- a Gaussian patch of width 0.2 on the
fault bent from :math:`x = 1/2`, its slip rising smoothly over 0.8 s from the
start -- and one time step for all meshes (0.5 ms), the four receivers of the
meshes of :math:`6^3`, :math:`8^3`, :math:`12^3` and :math:`16^3` differ from
the ones of :math:`24^3` by 1.6e-1, 5.7e-2, 3.1e-2 and 8.0e-3, and the slip
rates on the fault by 3.1e-2, 1.5e-2, 5.7e-3 and 2.1e-3, after 0.25 s -- as
much as with the same slip on the plane fault, since a receiver on a fault
takes the slip rate of the quadrature point nearest to it. Where the slip is
done within the run, the potency of the fault is the integral of the slip over
the bent plane to :math:`2.6 \cdot 10^{-6}`, on every mesh alike -- what is
left is the time integration of the slip rate.

What ``davschneller/config`` holds
-----------------------------------

A prototype of the per-cell-matrix formulation: a generated ``bootstrap`` kernel
that assembles a cell's ``M``, ``k``, ``kT``, ``r`` from quadrature weights times
a point-sampled metric, and the collocation data it needs. With the reference
mass matrix kept and the metric applied at the operator points, none of it is
needed here; it remains the reference for that formulation, should its flop
count win on some machine.
