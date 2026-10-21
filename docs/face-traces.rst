..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Face traces
===========

The corrector of a cell works from what its faces carry: after its predictor, a
cell writes the trace of its time integral on each of its faces, and the
corrector of a face reads that trace and the one the neighbour wrote on the same
face, applies the flux to the two and lifts the result into the cell. It reads no
field of the volume of either cell. This page says what a trace is, how the two
cells of a face agree on one, what local time stepping, dynamic rupture and the
device make of it, and what it was measured at.

The trace
---------

A face has a modal basis of its own, of the degree of the cell: 10 modes at
order 4, 21 at order 6, 36 at order 8. The trace of a field on face :math:`i` is
its restriction there, :math:`r_T(i)\,I` in the basis of the face (``rT``), exact,
since the face basis is orthogonal and ordered by degree. What the corrector of
the face does with the field of the cell, it can do with that.

The two cells sharing a face parametrise it differently: each sees the
parametrisation of the other reflected, by :math:`R = M_2^{-1} f_P(0)`, which is
its own inverse; the canonical vertex order makes the face orientation index
zero on every interior face. So every face has one parametrisation that both of
its traces are written in, the one of the cell the setup calls side 0 of the
face: of the two pairs of global id and face number, the one with the smaller
(``traceSideOf``), and on a fault the plus side. The other side writes its trace
reflected, :math:`R\,r_T(i)\,I` (``faceTraceMatrix``). Both cells tell from what
the mesh gives them about the face which side they are, without asking each
other, across ranks as well.

The corrector
-------------

Where a face carries its flux as a handful of scalars at its points -- in a
build with ``MATERIAL_NODAL`` for one of the isotropic media, elastic,
viscoelastic or acoustic -- the corrector forms the jump of the two traces, the
neighbour's less the cell's own, reads it at the points of the face
(``traceToFace``), turns it into the coordinates of the face, applies the
scalars there, turns the result back, integrates it against the face basis with
the weights of the rule the points are of, and lifts it into the cell. One set
of scalars does for both traces: in the strong form the corrector takes, the
part of the Godunov state of the cell's own trace is the part of the
neighbour's negated (``correctorScalars``). A face without a neighbour -- a
boundary or a fault -- reads a trace of zeros as the neighbour's
(``LTS::ZeroFaceTrace``) and takes its own part negated. A boundary whose ghost
state is a map of the cell's own state -- Dirichlet, and the free surface with
gravity -- folds that map into the matrix of the cell's own side, which no scalar
of the face carries: there the setup keeps the folded part alone in that matrix,
and the local integral lifts it from the time integral of the cell
(``appliesFoldApart``, ``foldedBoundaryFlux``).

That part takes the material of the cell, as the offset of the condition does,
while the cell's own part takes the material at the points of the face. Where
the material varies along such a face, the two do not cancel exactly for the
identity map, and a map that mirrors the tractions is not quite the free
surface: for an acoustic medium whose bulk modulus varies by a factor of 4 along
the surface and by up to 2 across a cell, the energies of the two differ by
1.7e-4, where those of the builds with and without ``MATERIAL_NODAL`` differ by
1.4e-2. Where the material is constant along the face, they are the same to
round-off.

Where the flux is one matrix per side -- in every build without
``MATERIAL_NODAL``, and for the poroelastic and the anisotropic medium in every
build -- each trace takes its own, :math:`t_O A^+ + t_N A^-`, and the lift goes
through ``traceLift``: :math:`\mathrm{rDivM}(i)\,M_2` for side 0 and
:math:`\mathrm{rDivM}(i)\,f_P(0)` for side 1, which takes the reflected trace
back. A face without a neighbour applies the matrix of the side outside the cell
to the trace of zeros; a free surface has no state there, and its matrix of that
side is zero (``hasExteriorOperator``).

The local integral writes the traces of a cell, the neighbour integral applies
the correctors of all of its faces.

Local time stepping
-------------------

What a cell writes on a face for the neighbour across it follows from how their
time steps relate (``LtsSetup::faceBuffer``):

- the same step: its trace of the step, which is its own trace;
- the neighbour takes a larger step: the traces of its steps, summed until the
  step of the neighbour is complete;
- the neighbour takes a smaller step: the traces of its time derivatives, which
  the neighbour integrates over its own interval. The :math:`k`-th derivative of
  the predictor of a cell of degree :math:`N` has the degree :math:`N-k` where
  the chain narrows them, and its trace the modes of that degree
  (``faceDerivatives``): over all derivatives as many modes as the cell has;
- a boundary face: nothing.

Across ranks, the cells exchange these per face: what a rank sends of a cell is
what the cells of the other rank read of it.

Dynamic rupture
---------------

The two sides of a fault face write the traces of their time derivatives, both
in the parametrisation of the plus side, and the fault evaluates them at its
points in time and in space with one matrix for both sides: on the device in one
kernel for all faces and both sides of a layer (``projectTracesToDR``). The fault
receivers and the static work of the energy output read the same traces.

Measured
--------

In a harness with the kernels of the production build (libxsmm and PSpaMM,
AVX2, one core), the two flux phases of a cell -- writing its traces and
applying its correctors -- against the volume form of the same flux, a local
flux from the time integral of the cell and a neighbour flux from the one of the
neighbour:

.. list-table::
  :header-rows: 1
  :widths: 40 15 15 15 15

  * -
    - order 4, DRAM
    - order 4, cache
    - order 6, DRAM
    - order 6, cache
  * - volume form
    - 4440 ns
    - 3740 ns
    - 19700 ns
    - 16900 ns
  * - traces and corrector
    - −27 %
    - −33 %
    - −40 %
    - −43 %

End to end, in SeisSol, with the kernels of a build that generates them through
Eigen (AVX2, one core): a planar wave through a constant elastic medium on 10368
cells, order 4, 627 steps, the time spent in the compute kernels, the median of
three runs each:

.. list-table::
  :header-rows: 1
  :widths: 40 20 20 20

  * -
    - volume form
    - traces and corrector
    -
  * - without ``MATERIAL_NODAL`` (one matrix per side)
    - 67.6 s
    - 61.4 s
    - −9 %
  * - ``MATERIAL_NODAL`` (scalars per point)
    - 128.4 s
    - 99.1 s
    - −23 %

The runs of each build spread by less than 2 % among themselves.

The traces and correctors give what the volume form gives, to round-off, on the
cubes of up to 1296 cells at order 4, on one rank and on two, with global and
with local time stepping: the energies to 1e-15 of their size for the elastic
medium in the default build and with ``MATERIAL_NODAL``, free surfaces and four
dynamic rupture runs included, whose fault receivers agree to 8e-14 and whose
fault output agrees to 1.1e-14; to 8e-16 for the acoustic medium, 4e-15 for the
poroelastic one and 1.8e-15 for the anisotropic one; and in the builds without
``MATERIAL_NODAL``, a Dirichlet face with an offset and the free surface with
gravity included, to 5.5e-16. Where the plain lift grows on a jump inside a
cell, a slab of curved cells and a viscoelastic one, the difference grows with
it, to 1.8e-14 and 4.3e-14 by :math:`t = 0.1`.

What is not done
----------------

- The map of a Dirichlet face and of the free surface with gravity at the points
  of the face. It takes the material of the cell (see above); a trace the map
  makes of the cell's own trace on such a face, which the scalars of the face
  would then take point by point as the trace of a neighbour, would take the
  material there, for a buffer per such face.
- The device beyond the compiler: the host side of its code compiles against
  the kernels the generator writes for CUDA (``sm_80``) for the elastic and the
  viscoelastic medium; nothing ran on a GPU. For the poroelastic medium the
  generator writes no kernels for a GPU.
