..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Material jumps inside a cell
============================

With ``MATERIAL_NODAL=ON`` a cell samples its material at the points its operator
is formed at, and a face at the points its flux is integrated at. The flux of a
face is the Riemann problem of the two materials there, times the coefficient
matrix of the cell's own, lifted into the cell. Where the material jumps inside
a cell, that is not stable: the lift carries the material of the face into the
whole cell, and the cell multiplies what reaches the other side of the jump with
the material there. A cell whose material spreads too far is therefore
*weight-adjusted* (Chan, Hewett and Warburton 2017): its faces lift the flux
without the material, and the cell applies the material once all of them are in,
at its sample points. This page says what that is, which cells take it, what it
costs, and what is not done yet.

When the plain lift grows
-------------------------

Measured in an acoustic prototype of the scheme -- the straight periodic mesh of
three by one by one cubes, eighteen cells, a layer whose edges cut cells, its
material jumping by a factor :math:`K` in the bulk modulus or in the density --
as the largest real part of the semi-discrete operator, which is round-off for a
stable one:

.. list-table::
  :header-rows: 1
  :widths: 22 14 16 16 16 16

  * - jump
    - impedance ratio
    - order 3
    - order 4
    - order 5
    - order 6
  * - bulk modulus, 4
    - 2
    - stable
    - stable
    - stable
    - stable
  * - bulk modulus, 9
    - 3
    - stable
    - stable
    - stable
    - stable
  * - bulk modulus, 16
    - 4
    - stable
    - 1.7e-2
    - stable
    - stable
  * - bulk modulus, 25
    - 5
    - stable
    - 1.2e-1
    - stable
    - stable
  * - bulk modulus, 100
    - 10
    - stable
    - 6.5
    - 6.6
    - stable
  * - density, 49
    - 7
    - stable
    - 2.5e-1
    - stable
    - stable
  * - density, 100
    - 10
    - stable
    - 7.2e-1
    - stable
    - stable

A transition over a tenth of a cell instead of a jump stays stable at order 4
even at a contrast of 100. Where the plain lift grows depends on the order and on
where the jump lies against the points: with the layer moved against the cells,
it grows at order 3 as well, from an impedance ratio of 5 on (see `Which cells,
and their time step`_). So the measurements set a bound to stay under rather than
a law.

The scheme
----------

The coefficient matrix of every direction is the inverse of the weight of the
medium times one without material, :math:`A_d = W^{-1} B_d`; for an elastic
medium :math:`W^{-1} = \operatorname{diag}(C, 1/\rho)`, and :math:`B_d` takes the
rate of the strain as the symmetric part of the velocity gradient and the rate of
the momentum as the divergence of the stress. :math:`C` takes the strain with its
tensor components, :math:`2\mu` on the shear diagonal: so the strain turns with a
face as the stress does, which is how a face turns what it lifts back into the
cell. With the engineering shear strains the two would differ by a factor of two
on the shear components. An acoustic medium has
:math:`W^{-1} = \operatorname{diag}(\lambda, 1/\rho, 1/\rho, 1/\rho)`.

A poroelastic medium takes the strain and the increment of fluid content as the
rates :math:`B_d` gives, and its weight couples them: :math:`W^{-1}` is the
undrained stiffness :math:`\bar C` on the strain, :math:`-M\alpha` between the
strain and the pressure and :math:`M` on the pressure, and on the velocities of
solid and fluid the inverse of the density matrix
:math:`\begin{pmatrix}\bar\rho & \rho_f \\ \rho_f & m\end{pmatrix}` in every
direction. :math:`B_d` takes the rate of the strain from the velocity of the
solid and the one of the fluid content from the velocity of the fluid, the
momenta from the divergence of the stress and the gradient of the pressure; the
drag between fluid and solid stays in the source term. The weight is symmetric
and positive definite, the energy of the medium.

An anisotropic medium has :math:`W^{-1} = \operatorname{diag}(C, 1/\rho)` and the
:math:`B_d` of the elastic one as well, :math:`C` its stiffness on the strain
with its tensor components. Its Voigt constants take a shear strain with the
engineering component, twice the tensor one, so the column of a shear strain
carries its constant twice: :math:`W^{-1}` is the symmetric Voigt matrix times
:math:`\operatorname{diag}(1, 1, 1, 2, 2, 2)` on the strain. That is not
symmetric where a normal strain reaches a shear stress, but self-adjoint
against the inner product of the strain energy, :math:`\sigma : \varepsilon`,
which counts every shear component twice as well, and positive definite there.
An isotropic medium couples no normal strain to a shear stress, and its weight
is symmetric either way. The kernels take the entries of :math:`W^{-T}`, as they
take the transposed coefficient matrices.

A medium with memory variables takes the flux and the weight of its base
material. Its memory variables take the strain rate times their relaxation
frequency, which carries no material, so the weight leaves the rows a face lifts
for them alone.

A medium states :math:`W^{-1}` in its own coefficients (``WeightEntry``) and the
coefficients of :math:`B_d` (``getUnweightedCoefficients``); the code generator
writes both into ``coefficients.h``.

A weight-adjusted cell

* lifts the flux of every face with :math:`B` instead of :math:`A`: the Godunov
  state of the two materials at a point of the face is the same, it is the matrix
  it is contracted with that changes (``fluxScalarsOfNode``). The scalars a face
  carries keep their positions. Two of them become zero: without material the
  flux of a face reaches no transverse normal strain, which only the weight turns
  into stress, inside the cell. The same goes for a face whose flux is one matrix
  per side, as the one of a poroelastic medium is (``liftedCoefficientMatrix``),
  for the lift of a fault face and for the boundary conditions that apply a
  matrix per face;
* collects those lifts in a buffer of its own (``LTS::FluxBuffer``), which the
  local integral clears and the neighbour integral finishes. The buffer has the
  shape of what the faces lift into: the field, or for the solver that keeps the
  memory variables apart (``SOLVER=linearckanelastic``) its extended field of the
  base quantities and a strain rate, which the local and the neighbour integral
  hand on to the field and the memory variables alike. So a weight-adjusted cell
  of that solver hands on its volume term in the local integral, and all of its
  fluxes, weighted, in the neighbour integral;
* adds at the end of the neighbour integral, before the plasticity, the buffer
  read at the sample points, times :math:`W^{-1}` there, projected back
  (``kernel::weightAdjust``); the rows beyond the quantities of the medium come
  in as they are;
* keeps its volume term as it is. On a straight-sided cell the volume term
  without material is a polynomial of the cell's degree, so applying the weight
  to it at the points is what the cell does already. A curved cell is
  weight-adjusted whatever its material does, weighs with its Jacobian
  determinant as well, and forms its volume term through the weight, in a split
  form (see :doc:`curvilinear-elements`).

A cell that is not weight-adjusted writes its fluxes into its field as before and
leaves the buffer alone; both kinds share their clusters. The Rusanov penalty
acts on the state and would have to be weighed with the material at the face, so
a weight-adjusted cell takes the Godunov flux on all of its faces.

On the device, the weight-adjusted cells of a layer form a batch of their own
(``KernelNames::WeightAdjust``): the local integral clears their buffers, every
kernel that lifts a face writes where the batch says the faces of a cell lift
into (``inner_keys::Wp::Id::FluxTarget``), and the neighbour integral applies the
weight to the batch (``gpu_weightAdjust``).

Which cells, and their time step
--------------------------------

``weightadjustedratio`` in ``&equations`` decides first: a cell may be
weight-adjusted where its impedances spread by more than that ratio over the
points its material is asked at -- the largest ratio between two of its P
impedances, or between two of its S impedances at points that carry shear
(``ImpedanceSpread``). An anisotropic medium gives the impedances of its
isotropic average there, from the means of its normal and of its shear
constants. It is 2 by default in a build that can do it, and 0 adjusts no cell.
The decision is taken where the time step of a cell is, and carried with it
through the partitioning.

A material that varies smoothly inside a cell keeps the plain lift stable where a
jump of the same spread does not. So such a cell is weight-adjusted only where
its impedances are not close to a smooth material (``ImpedanceJump``): the
logarithms of its P and of its S impedances at the samples are projected onto
the polynomials of one degree below the cell's, with the weights of the rule the
samples are of, and the cell is weight-adjusted where what the projection leaves
-- at the samples, and where the flux reads the material at the points of the
faces, at those points as well -- still spreads by more than
``weightadjustedjumpratio``, 2 by default. A jump between a face and the samples
nearest to it is at the points of that face alone. The projection overshoots at a
jump, so that spread can be larger than the one of the impedances themselves.
Where some points carry shear and others do not, the shear modulus itself jumps,
and the spread is unbounded. At order 2 the projection is onto the constants, and
what it leaves spreads as the impedances do. A ``weightadjustedjumpratio`` of 1
or less weight-adjusts every cell whose impedances spread by more than
``weightadjustedratio``, and with a ``weightadjustedratio`` below one then every
cell.

Measured in the prototype, on its mesh of eighteen cells: the layer with edges of
width 0.005, 0.015, 0.03 and 0.1 of the mesh -- a jump up to a third of a cell --
at three positions against the cells, its contrast 9 or 100 in the bulk modulus
or 25 in the density; and a smooth layer whose bulk modulus or density is
:math:`e^{A \sin 2\pi x}`, for :math:`A` up to 2 a spread of up to 5 inside a
cell. At orders 3 and 5 a part of these: the contrasts 25 and 100 in the bulk
modulus, the edges up to 0.03 wide, the smooth layer in the bulk modulus with
:math:`A` of 1 and 2. The operator with the cells a rule picks weight-adjusted,
the plain lift in the others: how many of the layers stay stable, how many cells
the rule picks, and how many of those in the smooth layers:

.. list-table::
  :header-rows: 1
  :widths: 40 20 20 20

  * - weight-adjusted where the impedances
    - order 3
    - order 4
    - order 5
  * - spread by more than 2
    - 20 of 20, 190, 18
    - 44 of 44, 414, 76
    - 20 of 20, 218, 24
  * - and the projection of one degree below leaves more than 2, at the samples
      and at the points of the faces
    - 20 of 20, 166, 0
    - 44 of 44, 222, 0
    - 20 of 20, 162, 0
  * - and the same projection leaves more than 1.5 at the samples alone
    - 20 of 20, 166, 2
    - 44 of 44, 228, 0
    - 20 of 20, 170, 0
  * - and the projection of the cell's own degree leaves more than 1.5 at the
      samples
    -
    - 41 of 44, 50, 0
    -

The plain lift grows in two of the layers at order 3 (edges 0.03 wide, contrast
25 and 100), five at order 4 (contrast 100, edges at most 0.015 wide) and one at
order 5. At the samples and at the points of the faces, the projection leaves at
most 1.16 of a smooth layer that spreads by 5 at order 4, 1.03 at order 5, and
1.81 at order 3, where it is the one of degree one, of the layer that spreads by
4.3 there; with 3 instead of 2, the plain lift grows in two of the layers at
order 3. At the samples alone, the projection leaves 1.53 of that layer at order
3, and two of its cells stay weight-adjusted with 1.5. The projection of the
cell's own degree follows the 24 samples of a jump at order 4 too closely: in
three of the five layers the plain lift grows in, the cells it leaves out make it
grow, with 1.25 instead of 1.5 as well.

A jump next to a face misses the samples: the same layer, an edge of it 0.003,
0.006 or 0.012 of the mesh inside the cells on either side of a face of them,
0.002 or 0.005 wide, contrast 9 or 100 in the bulk modulus or 25 in the density,
36 layers at orders 3 and 4; at order 5, where the samples come closer to the
faces, 0.001, 0.002 and 0.004 inside, 0.002 wide, contrast 100 in the bulk
modulus or 25 in the density, 12 layers. The plain lift grows in one of them at
order 3, eleven at order 4 and none at order 5. The projection at the samples
alone, with 1.5, keeps all but one stable: at order 4, an edge 0.006 inside a
cell, 0.002 wide, contrast 100, between a face and the samples nearest to it,
where it leaves out two cells the spread takes. At the samples and at the points
of the faces, with 2, every one of them stays stable, with 216, 300 and 120
cells weight-adjusted (218, 300 and 120 by the spread alone); with 3 instead of
2 the plain lift grows in one of them at order 4.

A weight-adjusted cell lifts the dissipation of the Riemann problem of the face
and multiplies it with the stiffest material inside the cell, and pays with its
time step. Measured in the same prototype as the spectral radius over the one of
the exact weighted mass matrix, which has no penalty:

.. list-table::
  :header-rows: 1
  :widths: 28 18 18 18 18

  * - jump, impedance ratio
    - order 3
    - order 4
    - order 5
    - order 6
  * - bulk modulus, 2
    - 1.02
    - 1.21
    - 1.01
    - 1.01
  * - bulk modulus, 5
    - 1.08
    - 2.63
    - 1.43
    - 1.28
  * - bulk modulus, 10
    - 1.33
    - 4.83
    - 2.32
    - 2.12
  * - density, 5
    - 1.25
    - 1.85
    - 1.07
    - 0.96
  * - density, 10
    - 2.29
    - 3.45
    - 1.80
    - 1.42

A transition over a tenth of a cell costs at most 1.4 at a contrast of 100.
:math:`(1 + r)/2` bounds all of these, but in three dimensions a cell needs more.
On the periodic cube of the next section (CFL 0.5), a slab whose jumps cut the
cells stays stable up to a time step of :math:`c/r` times the one its wave
speeds permit, :math:`c` the same for every :math:`r` but not for every order.
The largest :math:`c` that ran through to :math:`t = 1`, and the smallest that
grew (a forced time step factor :math:`c/r`):

.. list-table::
  :header-rows: 1
  :widths: 34 8 14 14 14 14

  * - jump
    - :math:`r`
    - order 3
    - order 4
    - order 5
    - order 6
  * - acoustic, bulk modulus by 9
    - 3
    - 3.0 at factor 1
    - 1.5 / 1.8
    - 1.8 / 2.0
    - 1.5 / 1.65
  * - acoustic, bulk modulus by 25
    - 5
    - --
    - 1.65 / 2.0
    - --
    - --
  * - acoustic, bulk modulus by 100
    - 10
    - 3.0 / 3.5
    - 1.7 / 1.8
    - 2.0 / 2.1
    - 1.63 / 1.76
  * - acoustic, bulk modulus by 400
    - 20
    - 3.0 / 3.4
    - 1.6 / 1.9
    - 2.0 / 2.2
    - --
  * - elastic, :math:`\lambda` and :math:`\mu` by 100
    - 10
    - --
    - 2.0 / 2.2
    - --
    - --

The time step of a weight-adjusted cell is shortened by :math:`c/r`, with a
margin of 1.5 on the smallest :math:`c` that ran through at its order: 2 at
order 3, 1 at orders 4 and 6, 1.2 at order 5. Orders 2, 7 and 8 are not measured
and take the smallest, 1 (``weightAdjustedTimestepScale``). The shear waves of an
elastic medium count the same way, weighed with their speed against the fastest
P wave. Where the material only varies smoothly, that is more than needed.

End to end
----------

Measured in SeisSol: the periodic cube :math:`[0, 2]^3` of six by six by six
cubes, six cells each, the vertices inside moved (1296 straight-sided cells);
elastic and acoustic, order 4, local time stepping at rate 2, CFL 0.5; a planar
wave through a slab :math:`0.6 < x < 1.4` whose faces cut cells, the density one
everywhere. The runs of this section choose the cells by the spread of their
impedances alone, ``weightadjustedjumpratio = 1``, and weight-adjust every cell
with ``weightadjustedratio = 0.5`` on top of it.
The strain energy the energy output reports is no measure here: it interpolates
the material between the samples, which overshoots at a jump inside a cell, and
is negative from the start. With one density the kinetic energy is exact:

.. list-table::
  :header-rows: 1
  :widths: 40 20 20 20

  * - slab
    - weight-adjusted cells, their time step
    - kinetic energy at :math:`t = 1`, plain
    - weight-adjusted, over :math:`0 \le t \le 1`
  * - :math:`\lambda` and :math:`\mu` 100 times as large: P and S impedances by 10
    - 414 of 1296, 0.1
    - 1.8e45, from 1.05
    - 0.42 to 1.05
  * - the same, its faces at an angle to the cells (normal :math:`(0.6, 0.8, 0)`)
    - 552 of 1296, 0.1
    - 5.0e39, from 1.05
    - 0.40 to 1.05
  * - :math:`\mu` 100 times as large, :math:`\lambda + 2\mu` the same: S impedances by 10
    - 414 of 1296, 0.14
    - 4.2e67, from 0.75
    - 0.30 to 0.75
  * - acoustic, the bulk modulus 100 times as large: impedances by 10
    - 414 of 1296, 0.1
    - 5.2e23, from 0.40
    - 0.24 to 0.40

A slab a hundred times as dense, its impedances by 10 as well, grows with neither
by :math:`t = 1`. Where the material varies smoothly -- a cosine bump by 2 in the
density and by 4 in :math:`\lambda` and :math:`\mu` over :math:`0.5 < x < 1.5` --
the weight-adjusted cells are as accurate as the plain ones: with all of them
weight-adjusted, eight receivers differ from a run on cells half as large by 2.25e-2
in the stresses (relative, over the whole run), and by 2.21e-2 with none.

With the default, ``weightadjustedjumpratio = 2``, the first slab weight-adjusts
the same 414 cells, and its energies are those of the run with 1 bit for bit. A
material that varies smoothly by as much -- :math:`\lambda` and :math:`\mu` by
:math:`\exp(3 \sin 2 \pi x)` over :math:`0 < x < 1`, the impedances by up to 15
across a cell and by 20 over the layer -- has 648 cells whose impedances spread
by more than 2. With 1 all of them are weight-adjusted, at 0.079 times the time
step their wave speeds permit; with the default none is, and the run takes 19 s
of compute instead of 74 s. Its kinetic energy stays between 0.76 and 1.74 to
:math:`t = 1`, that of the weight-adjusted run between 0.75 and 1.74, and the two
differ by 1.2 % of it at most.

The shorter time step is needed: at the one the wave speeds permit, the
weight-adjusted cells of the first and the third slab reach a kinetic energy of
1e183 and 1e96 by :math:`t = 0.05`; the limits between are in the previous
section.

With every cell weight-adjusted on a constant material, four by four by four
cubes, the errors against the planar wave are the plain run's to 4e-15, on one
rank and on two, and to 1e-15 for an acoustic medium; with the default ratio no
cell is, no buffer is allocated, and the run is the plain one bit for bit. The
first slab, started from a pulse given by easi instead -- the planar wave takes
its material from the first cell of each rank -- gives the same kinetic energy on
one rank and on two to 9e-16.

On curved cells -- the same cube, its cells of order two, ``CURVILINEAR`` with
the assembled operator -- every cell, being curved, is weight-adjusted with its
Jacobian determinant and forms its volume term in the split form (see
:doc:`curvilinear-elements`), and the 406 cells the slab cuts take the time
step of their impedances as well: the first slab stays between 0.40 and 1.05 to
:math:`t = 1` and ends at 0.50.

A viscoelastic medium with three mechanisms, :math:`Q_P = 40` and :math:`Q_S = 20`
everywhere at a central frequency of 1 and a ratio of 100, with the solver that
keeps the memory variables apart (``SOLVER=linearckanelastic``): the first slab
grows from 0.089 to 6.5e30 by :math:`t = 0.5` with the plain lift; with its 414
cells weight-adjusted it stays at or below 0.089 to :math:`t = 1`, decaying to
0.028. With every cell weight-adjusted on a constant material, the energies --
those of the memory variables among them -- are the plain run's to 4e-16, and so
are the errors against the planar wave in the digits the log gives.

A viscoacoustic medium with three mechanisms and :math:`Q_P = 40` everywhere, the
bulk modulus 100 times as large inside the slab, with the same solver: it grows
from 0.048 to 4.5e11 by :math:`t = 0.5` with the plain lift, and stays at or below
0.048 to :math:`t = 0.5` with its 414 cells weight-adjusted. With every cell
weight-adjusted on a constant material, the energies are the plain run's to
5e-16.

A poroelastic medium, the slab with :math:`\lambda`, :math:`\mu` and the bulk
modulus of the solid 100 times as large, the fluid, the porosity and the
densities the same everywhere: it grows from 1.26 to 2.3e40 by :math:`t = 0.25`
with the plain lift; with its 410 cells weight-adjusted, at 0.13 times the time
step their wave speeds permit, it stays between 0.58 and 1.26 to :math:`t = 1`.
With every cell weight-adjusted on a constant material, the energies are the
plain run's to 5e-16, and so are the errors against the planar wave in the
digits the log gives.

An anisotropic medium -- a stiffness that couples a normal strain to a shear
stress (:math:`c_{16}`) and two shear strains (:math:`c_{45}`), every constant
100 times as large inside the slab, the density one everywhere: it grows from
6.0 to 7e180 in kinetic energy by :math:`t = 1` with the plain lift, and stays
between 2.9 and 6.0 to :math:`t = 1` with the 410 cells weight-adjusted that the
impedances of its isotropic average pick, at 0.1 times the time step. With every
cell weight-adjusted on the medium outside the slab everywhere, the energies are
the plain run's to 3e-16, with a free surface as well, on one rank and on two.

Cost
----

The kernel that applies the weight, per weight-adjusted cell and step, against
the predictor, the volume term and the eight fluxes of the cell (elastic, double
precision, ``hsw``, thousands of nonzero / hardware operations):

.. list-table::
  :header-rows: 1
  :widths: 20 26 26 28

  * - order
    - cell and step
    - weight
    - share
  * - 4
    - 118.1 / 136.7
    - 15.7 / 17.9
    - 13 %
  * - 6
    - 792.0 / 909.0
    - 145.8 / 161.2
    - 18 %
  * - 8
    - 3666.1 / 4017.2
    - 726.1 / 786.3
    - 20 %

An acoustic cell at order 4: 6.8 / 7.8 on 46.3 / 51.8. A viscoelastic one with
three mechanisms: 16.1 / 18.2 on 323.9 / 393.3 with the memory variables kept
apart (5 %), 16.6 / 18.5 on 593.7 / 680.7 with them in the field (3 %). A
poroelastic one: 23.3 / 26.5 on 690.1 / 1011.5 (3 %), most of which is its
space-time predictor.

Timed on one core (a Xeon at 2.1 GHz, the kernels for ``hsw``) over 4096 cells
with operands of their own, the weight takes 0.9 µs against 11.5 µs for the rest
of an elastic cell's step at order 4 (8 %), and 9.4 µs against 79 µs at order 6
(12 %): it is a product with dense matrices, which runs faster than the
predictor.

What a run pays is mostly the time step. In the first slab above, where a third of
the cells is weight-adjusted and steps at a tenth of what its wave speeds permit,
the clusters update 6.9 times as many cells per time as in the plain run, and the
run took 6.4 times as long; with the S impedances alone by 10, 3.9 times as many,
and for the smooth bump with every cell weight-adjusted 1.6 times as many. The
viscoelastic and the viscoacoustic slab update 6.9 times as many cells per time,
the poroelastic one 5.2 times as many.

Every cell of a rank carries the buffer where any of them is weight-adjusted, as
large as the field its faces lift into: 1.4 kB for an elastic cell at order 4,
4.0 kB at order 6, 2.4 kB for a viscoelastic one at order 4 with the memory
variables kept apart, in double precision. A rank without such a cell has none.

What is not done
----------------

- A jump in the anisotropy itself. An anisotropic cell is chosen by the
  impedances of its isotropic average, which a turn of the axes of the medium
  changes little or not at all, and its time step is shortened by the scale
  measured for isotropic media, with the spread of that average.
- The memory variables with ``SOLVER=linearck``. Its kernels take the weight as
  the others do and compile, but a viscoelastic build with that solver does not,
  weight or not (the energy output has no viscoelastic part for it), so none of
  it ran.
- The device beyond the compiler. The host side of its code compiles against the
  kernels generated for CUDA (``sm_80``), and the weight kernels the generator
  writes compile to machine code for ``sm_80`` (clang to PTX, ``ptxas``); nothing
  ran on a GPU. The generator for AMD GPUs stops in TensorForge for the nodal
  elastic build, also without the weight kernel.
- The time step at orders 2, 7 and 8, which take the smallest :math:`c` of the
  orders measured.
- ``weightadjustedjumpratio`` at orders 2 and 6 to 8: its default of 2 is
  measured at orders 3 to 5.

Checks
------

Unit checks: the operator of a medium is its inverse weight times the one without
material, in every direction, for solids, a fluid in an elastic build, acoustic,
anisotropic and poroelastic media, exactly, in Python on the declarations and in
C++ against the matrices the materials write themselves, and for the media with
memory variables against the operator of their solver, the memory variables
included; every weight is self-adjoint against the inner product of the strain
energy, and the isotropic ones symmetric; the weight of a poroelastic and of an
anisotropic medium is positive definite there; the flux without
material decomposes into the scalars of a face as the one with the material
does, between solids, fluids and both, and at a free surface; the kernel that
applies the weight does what its matrices say where the material differs from
sample to sample, a fluid among the samples, and leaves the rows beyond the
quantities of the medium as they are; with one material for the cell, what a
weight-adjusted cell adds -- the volume term, the correctors of a regular, a
free-surface and a fault face, the lift of a fault face, a face
that carries a matrix per side -- is what a cell adds that lifts with the
material, for a solid and for a fluid, into the field and into the extended field
of the solver that keeps the memory variables apart; the spread of the
impedances and the time step it gives at every order; and what the projection of
a cell leaves of its impedances: nothing where their logarithms are a polynomial
of degree one, from order 3 on, or two, from order 4 on, more than a factor of
five of a jump by 10 across a plane through the cell, and an unbounded ratio
where the shear modulus vanishes at some samples only.

The device code: its host side compiles against the kernels the generator writes
for CUDA, for the elastic, the acoustic and the viscoelastic media with both
solvers, and the weight kernels compile to machine code for ``sm_80``.

End to end: the runs above -- a weight-adjusted run on a constant material
against the plain one, elastic, acoustic, viscoelastic, viscoacoustic and
poroelastic; the slabs against the plain lift and against longer time steps, at
orders 3 to 6, and the first one on curved cells; a smooth material against
cells half as large; one rank against two.
