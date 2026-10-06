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

Which cells are curved, the mesh decides. A reader hands them over through
``MeshReader::setCurvedGeometry``: the nodes of every local cell on the
equispaced lattice of the geometry order
(``IsoparametricTransform::latticeNodes``; order two adds the six edge midpoints
to the vertices), the vertices first and in the vertex order the reader keeps.
The call checks that the vertices are the cell's and that the cell is not
turned inside out anywhere on a lattice of twice its order. **No reader does
this yet**; a higher-order mesh through PUMgen and the PUML format is the next
step, in the mesh toolchain. A build without ``CURVILINEAR`` refuses a curved
mesh.

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
    - would be per point
    - The fault rotates its state with one rotation per face and lifts it with
      one scale. Not done; refused on a curved mesh.
  * - Boundary conditions with a ghost state
    - would be per node
    - Analytical, Dirichlet and gravity state a ghost cell and apply one matrix
      per face. Not done; refused on a curved mesh. Free surface and outflow are
      part of the flux and work.

What is not done
----------------

- A mesh reader that delivers curved cells (PUMgen, PUML), and with it the
  transforms for the time step at partitioning time (``CellToVertexArray``,
  which ``fromPUML`` builds straight-sided).
- Dynamic rupture and the boundary conditions with a ghost state on curved
  faces, see above.
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
face included.

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

What ``davschneller/config`` holds
-----------------------------------

A prototype of the per-cell-matrix formulation: a generated ``bootstrap`` kernel
that assembles a cell's ``M``, ``k``, ``kT``, ``r`` from quadrature weights times
a point-sampled metric, and the collocation data it needs. With the reference
mass matrix kept and the metric applied at the operator points, none of it is
needed here; it remains the reference for that formulation, should its flop
count win on some machine.
