..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Curvilinear elements
====================

Nothing in the code does this yet. This page records what a curved element
asks of the discretisation, which of those pieces already exist, and the three
things about it that are easy to get wrong. It is written so that whoever picks
the work up does not have to rediscover the measurements.

File references are to ``davschneller/nodal`` at ``a1cc532e5``; line numbers
are a hint, the symbol is the thing to look for.

What changes
------------

A cell is given by a map :math:`x(\xi)` from the reference tetrahedron. For an
affine map the Jacobian :math:`J = \partial x / \partial \xi` is one matrix for
the whole cell; for a curved map it is a function of :math:`\xi`. Everything
else follows from that.

The volume term of the DG formulation is

.. math::

  \int_T \frac{\partial \varphi_k}{\partial \xi_e}\,
         G_{ed}(\xi)\, A_d\, q \; \mathrm{d}\xi ,
  \qquad
  G = \det(J)\, J^{-1} ,

and the mass matrix is :math:`M_{kl} = \int_T \varphi_k \varphi_l \det J \,
\mathrm{d}\xi`. With an affine map both :math:`G` and :math:`\det J` are
constants and pull out of the integral, which is why the reference matrices in
``codegen/matrices/aderdg-N.xml`` suffice and a cell carries only the three
rows of :math:`J^{-1}`. With a curved map they do not.

For an isoparametric P2 map the quantities involved are polynomials: :math:`J`
is linear, :math:`\det J` is cubic, and :math:`G`, being the cofactor matrix,
is quadratic. Nothing is approximated by writing them down; what changes is
that they can no longer be factored out.

What is already in place
------------------------

Considerably more than one would expect, because three unrelated pieces of work
left exactly the right seams.

.. list-table::
  :header-rows: 1
  :widths: 24 30 46

  * - Piece
    - Where
    - What it gives
  * - Geometry abstraction
    - ``src/Geometry/CellTransform.h``, ``src/Geometry/FaceTransform.h``
    - ``CellTransform::refToSpaceJacobian(point)`` is virtual and takes a
      reference coordinate; ``FaceTransform::normal(faceCoord)`` is
      point-dependent and unnormalised, with its norm being the surface
      Jacobian. A curved cell is a new subclass, not a new interface.
  * - The metric as an operand
    - ``referenceGradients(dim)`` in ``codegen/kernels/aderdg/aderdg.py``
      (~335), ``LocalIntegrationData::referenceGradients`` in
      ``src/Initializer/Typedefs.h`` (~61)
    - A cell carries the three rows of :math:`J^{-1}` apart from the material,
      instead of a star matrix with both folded together. Filled in
      ``src/Initializer/Model/CellLocalMatrices.cpp`` (~299) from
      ``refToSpaceJacobianInverse(ReferenceBarycenter)``, with a comment
      marking that an affine map is what permits a single point there.
  * - An operator formed per sample point
    - ``MATERIAL_OPERATOR=assembled`` → ``starAtPoint(dim)``,
      ``aderdg.py`` (~456, ~851)
    - ``starAtPoint[dim][n,q,p]`` is an operator per point. A point-varying
      metric needs exactly this shape; the factored form's economy is that the
      Jacobian rows fold once per cell, which a curved cell denies.
  * - A time integration that tolerates a degree-raising operator
    - ``codegen/kernels/aderdg/linearck.py`` (~309),
      ``codegen/kernels/aderdg/stp.py``
    - The derivative chain keeps every mode instead of narrowing, and the
      space-time predictor is a Picard fixed point instead of a block sweep over
      degrees. See `The degree cascade does not survive`_ for why this is not
      optional.
  * - A flux that varies along a face
    - ``fluxScalarsOfNode`` in ``CellLocalMatrices.cpp``,
      ``nodalFlux`` in ``aderdg.py`` (~647)
    - The flux operator of a face is carried as scalars per face node, and the
      material of both sides is read per node. The surface Jacobian already
      rides on those scalars through ``fluxScale``.
  * - Point location for output
    - ``CellTransform::spaceToRef``
    - Receivers and output points need :math:`\xi` from :math:`x`; for a curved
      map this is a Newton iteration, and the virtual method is the place for it.

Three things that are easy to miss
----------------------------------

All three are measured, not argued. The scripts are small and self-contained;
they use SeisSol's own basis and quadrature, read from the collocation data on
``davschneller/config``.

The degree cascade does not survive
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ADER time integration rests on the spatial operator lowering the polynomial
degree: the derivative of a degree-:math:`n` field is degree :math:`n-1`, so the
derivative chain narrows and the space-time predictor can be solved one degree
block at a time. A polynomial metric destroys that property, because
:math:`\partial \varphi_k / \partial \xi_e \, G_{ed}` raises the degree by as
much as :math:`G` carries.

Measured on ``kDivM``: the largest entry in the degree blocks that a
degree-lowering operator must leave empty, relative to the largest entry
overall.

.. list-table::
  :header-rows: 1
  :widths: 50 50

  * - Cell
    - Outside the degree-lowering structure
  * - affine, :math:`J = I`
    - 1.9e-16
  * - affine, general :math:`J`
    - 2.2e-16
  * - P2-curved, :math:`\varepsilon = 0.05`
    - 7.7e-03
  * - P2-curved, :math:`\varepsilon = 0.2`
    - 4.2e-02

So a curved build has to take the non-narrowing derivative chain and the
fixed-point predictor. Those exist, but they are switched on by
``MATERIAL_NODAL`` — by a statement about the *material*. Curvilinear is a
second, independent reason for the same switch, and the condition wants to be
named after what it asserts ("the operator does not lower the degree") rather
than after one of its causes. This is the one code generator change beyond the
matrices themselves.

The mass matrix becomes dense
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The modal basis is orthogonal on the reference tetrahedron, so :math:`M` is
diagonal there. With :math:`\det J` under the integral it is not. Measured
occupancy, order 4: **5 % → 100 %**, at an essentially unchanged condition
number (84 → 80).

:math:`M^{-1}` is therefore a dense matrix per cell, and it is folded into more
places than the volume term. See `Where the inverse mass matrix is implicit`_.

The quadrature is one degree short
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The 44-point volume rule shipped for order 4 is exact to degree 8 — ample for
the affine integrands (:math:`\varphi \varphi` is degree 6). The curved mass
matrix integrand is :math:`\varphi \varphi \det J`, degree 9 for a P2 map, so it
is under-integrated. Measured against an exact rule:

.. list-table::
  :header-rows: 1
  :widths: 50 50

  * - Cell
    - Error in :math:`M`, relative
  * - affine
    - 2.7e-15
  * - P2-curved, :math:`\varepsilon = 0.05`
    - 3.2e-07
  * - P2-curved, :math:`\varepsilon = 0.2`
    - 2.0e-05

This is below the discretisation error, so it will not visibly break anything —
which is exactly why it is worth writing down. It costs exact conservation and
it is avoidable: the rule needs :math:`2p + \deg \det J` rather than
:math:`2p`, which for a P\ :sub:`g` map is :math:`2p + 3(g-1)`.

Two formulations, and the storage they cost
--------------------------------------------

There are two ways to carry a curved cell's operator, and they compute the same
thing. The choice is a memory-against-flops trade and deserves a measurement on
the target machine rather than a preference.

**Per-cell matrices.** Integrate the metric into the cell's own copies of the
six reference families in ``codegen/matrices/aderdg-N.xml`` — ``kDivM``,
``kDivMT``, ``rDivM``, ``fMrT``, ``rT``, ``fP`` — and let the cell's star be the
raw directional matrices :math:`A_d^T`. The kernels then look exactly as they do
today; only their operands become per-cell. This is what ``davschneller/config``
implements.

**Per-point metric.** Keep the reference matrices and carry the metric at the
sample points, applied the way the material coefficients already are, with one
dense :math:`M^{-1}` per cell for the projection back. This is the shape
``MATERIAL_OPERATOR=assembled`` already generates.

Storage per cell, f64, elastic:

.. list-table::
  :header-rows: 1
  :widths: 16 14 24 24 22

  * - Order
    - DOFs ``Q``
    - Per-cell matrices
    - :math:`M^{-1}` only
    - Metric at ``nb`` points
  * - 4
    - 1.4 kB
    - 40.6 kB
    - 3.1 kB
    - 1.4 kB
  * - 6
    - 3.9 kB
    - 271.0 kB
    - 24.5 kB
    - 3.9 kB

A factor of roughly nine at order 4 and ten at order 6, against more work per
application. Note also that the per-cell matrices are dense in the curved case —
the degree structure that makes the affine ``kDivM`` sparse is the same structure
the metric destroys.

The face rotation is a cost in either formulation. A curved face has a normal
that varies along it, so ``T`` and its inverse become per node:
``faceRotation[side]`` in ``src/Initializer/Typedefs.h`` (~90) and
``rotation["qk"]`` / ``rotation["pl"]`` in ``nodalFlux``. Stored as a matrix per
node that is 14.1 kB per cell at order 4 and 29.5 kB at order 6, against 1.4 kB
per face today. ``T`` has 45 stored entries but only three degrees of freedom
per node, so keeping the normal per node and building the rotation in the kernel
is about fifteen times cheaper in memory and correspondingly more work in the
kernel.

Where the inverse mass matrix is implicit
------------------------------------------

This is the part that is easy to underestimate. :math:`M^{-1}` is not a matrix
the code multiplies by; it is pre-multiplied into the checked-in reference
matrices, and those reach the kernels through the generated pool rather than
through an operand the host sets. Every one of these has to learn about a
per-cell mass matrix.

.. list-table::
  :header-rows: 1
  :widths: 26 32 42

  * - Matrix
    - Named in
    - Reaches
  * - ``kDivM``, ``kDivMT``
    - ``aderdg/{aderdg,linearck,linearckanelastic,stp}.py``
    - volume term, derivative chain, space-time predictor
  * - ``rDivM``, ``fMrT``
    - ``aderdg/{aderdg,linearck,linearckanelastic}.py``
    - local and neighbour flux
  * - ``project2nFaceTo3m``
    - ``nodalbc.py``, the three solvers
    - nodal boundary conditions, the nodal flux lift
  * - ``V3mTo2nTWDivM``
    - ``dynamic_rupture.py``
    - the dynamic rupture flux
  * - ``M3inv``
    - ``point.py``
    - point sources
  * - ``M2inv``
    - ``aderdg/aderdg.py``
    - the face reparametrisation check
  * - ``M2``, ``MV2nTo2m``
    - ``surface_displacement.py``
    - free-surface displacement — a face integral, so it wants the curved
      surface Jacobian rather than the volume mass matrix

Three more places integrate or evaluate over a cell without naming one of these,
and want checking rather than assuming: plasticity (``plasticity.py``, which
transforms modal to nodal with ``evalAtQP``/``vInv`` and accumulates plastic
strain), the energy output, and the initial-condition projection. The output
projection in ``vtkproject.py`` evaluates basis functions at points and needs no
mass matrix, but it does need ``spaceToRef`` for a curved map.

A route through it
------------------

Ordered so that each step can be checked against the affine result before the
next one starts. An affine transform evaluated at many points must reproduce,
bit for bit, what one evaluated at the barycenter produces — that comparison is
the test harness for the whole series.

1. **Bring the branch up to master.** ``davschneller/nodal`` has the geometry
   abstraction and the refactored boundary conditions through its merge base, but
   master keeps moving.

2. **Give the metric a point index.** ``referenceGradients`` gains a point
   index, ``MATERIAL_OPERATOR=assembled`` becomes the required form, and
   ``CellLocalMatrices.cpp`` fills it from ``refToSpaceJacobianInverse`` at each
   sample point instead of at the barycenter. Check: an ``AffineTransform``
   sampled at every point against the current path.

3. **Separate the degree condition from the material.** The non-narrowing
   derivative chain and the fixed-point predictor are selected by
   ``self.nodalMaterial`` today; curvilinear needs the same behaviour for a
   different reason. One predicate, two causes.

4. **A per-cell inverse mass matrix, and then the survey above.** Expect this to
   be the longest step — not the kernels, but the places :math:`M^{-1}` has been
   folded into.

5. **The face rotation per node.** ``FaceTransform::faceAlignedBasis()`` is the
   one accessor on that interface that does not yet take a point, while
   ``normal(input)`` and ``surfaceJacobian(input)`` do. Then
   ``MeshTools::normalAndTangents`` in ``CellLocalMatrices.cpp`` (~365) gives way
   to the transform, and ``fluxScale`` (~393) becomes a per-node quantity rather
   than a single ratio of face area to cell volume.

6. **A curved transform, and a mesh that carries one.** Only now does a
   non-affine ``CellTransform``/``FaceTransform`` subclass become useful, and
   with it the question of where curved geometry enters — PUMgen, the PUML
   format, and a higher-order mesh from the mesher.

Steps 2, 3 and 5 are each small and well-bounded. Step 4 is a survey. Step 6 is
a separate project in the mesh toolchain.

Two notes on the sample points. The ``nb`` set has one point per basis function
and its face traces are the two-dimensional nodal set, so the same samples serve
the volume metric and the face rotation; the ``ip`` set integrates products
further but has no point on a face, so a curved build wants ``nb`` unless the
face nodes are sampled separately. And a curved geometry raises the degree under
every integral, which is the quadrature point above.

What ``davschneller/config`` holds
-----------------------------------

The branch is a prototype of the per-cell-matrix formulation. Its merge base is
well behind master and its own ``CellTransform`` has been superseded by the one
in ``src/Geometry``, so it is more useful read than merged. What is in it and
nowhere else:

- ``codegen/kernels/elementwise.py`` — a generated ``bootstrap`` kernel that
  assembles a cell's ``M``, ``k``, ``kT``, ``r`` from quadrature weights times a
  point-sampled metric. This is the piece worth taking.
- ``codegen/matrices/elemwise-collocate-pN.json`` — basis values and
  derivatives at volume and face quadrature points, which is what any
  formulation needs to sample a metric. Regenerate these a degree higher; see
  the quadrature note above.
- ``MatrixBootstrap::sampleBasis`` — ``det(J)`` and ``det(J) J^{-1}`` per
  quadrature point, with a rescaling of the weights for conditioning that
  cancels later.
- ``Config::GlobalElementwise`` and the ``LTS`` variables that hold the per-cell
  families, as a worked example of the plumbing.

Open decisions
--------------

- Which formulation. The storage table above is the argument; the flop count on
  the target machine is the other half, and it has not been measured.
- Whether the face rotation is stored per node or rebuilt in the kernel from a
  stored normal.
- Which geometric order to support. P2 is the natural first step and fixes the
  polynomial degrees used throughout this page; anything higher scales the
  quadrature requirement as :math:`2p + 3(g-1)`.
- Whether a curved build must also be a ``MATERIAL_NODAL`` build. The two share
  the per-point operator machinery, but a curved cell with a constant material is
  a legitimate and cheaper configuration, and keeping it separate is what step 3
  above is for.
