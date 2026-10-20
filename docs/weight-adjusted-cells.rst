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
where the jump lies against the points, so the measurements set a bound to stay
under rather than a law.

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
on the shear components. A medium states :math:`W^{-1}` in its own coefficients
(``WeightEntry``) and the coefficients of :math:`B_d`
(``getUnweightedCoefficients``); the code generator writes both into
``coefficients.h``.

A weight-adjusted cell

* lifts the flux of every face with :math:`B` instead of :math:`A`: the Godunov
  state of the two materials at a point of the face is the same, it is the matrix
  it is contracted with that changes (``fluxScalarsOfNode``). The scalars a face
  carries keep their positions. Two of them become zero: without material the
  flux of a face reaches no transverse normal strain, which only the weight turns
  into stress, inside the cell. The same goes for the lift of a fault face and for
  the boundary conditions that apply a matrix per face;
* collects those lifts in a buffer of its own (``LTS::FluxBuffer``), which the
  local integral clears and the neighbour integral finishes;
* adds at the end of the neighbour integral, before the plasticity, the buffer
  read at the sample points, times :math:`W^{-1}` there, projected back
  (``kernel::weightAdjust``);
* keeps its volume term as it is. On a straight-sided cell the volume term
  without material is a polynomial of the cell's degree, so applying the weight
  to it at the points is what the cell does already. On a curved cell it is not;
  in the prototype the two then grow alike, by what a curved cell grows with one
  material as well.

A cell that is not weight-adjusted writes its fluxes into its field as before and
leaves the buffer alone; both kinds share their clusters. The Rusanov penalty
acts on the state and would have to be weighed with the material at the face, so
a weight-adjusted cell takes the Godunov flux on all of its faces.

Which cells, and their time step
--------------------------------

``weightadjustedratio`` in ``&equations`` decides: a cell is weight-adjusted where
its impedances spread by more than that ratio over the points its material is
asked at -- the largest ratio between two of its P impedances, or between two of
its S impedances at points that carry shear (``ImpedanceSpread``). It is 2 by
default in a build that can do it; 0 adjusts no cell, and a ratio below one every
cell. The decision is taken where the time step of a cell is, and carried with
it through the partitioning.

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
everywhere.
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

The shorter time step is needed: at the one the wave speeds permit, the
weight-adjusted cells of the first and the third slab reach a kinetic energy of
1e183 and 1e96 by :math:`t = 0.05`; the limits between are in the previous
section.

With every cell weight-adjusted (``weightadjustedratio = 0.5``) on a constant
material, four by four by four cubes, the errors against the planar wave are the
plain run's to 4e-15, on one rank and on two, and to 1e-15 for an acoustic
medium; with the default ratio no cell is,
no buffer is allocated, and the run is the plain one bit for bit. The first slab,
started from a pulse given by easi instead -- the planar wave takes its material
from the first cell of each rank -- gives the same kinetic energy on one rank and
on two to 9e-16.

On curved cells -- the same cube, its cells of order two, ``CURVILINEAR`` with
the assembled operator -- the first slab grows from 1.05 to 4.4e59 by
:math:`t = 1` with the plain lift and stays between 0.40 and 1.05 with 406 cells
weight-adjusted; with every cell weight-adjusted on a constant material, the
errors against the planar wave are the plain run's to 4e-15.

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

An acoustic cell at order 4: 6.8 / 7.8 on 46.3 / 51.8. Timed on one core (a Xeon at
2.1 GHz, the kernels for ``hsw``) over 4096 cells with operands of their own, the
weight takes 0.9 µs against 11.5 µs for the rest of the cell's step at order 4
(8 %), and 9.4 µs against 79 µs at order 6 (12 %): it is a product with dense
matrices, which runs faster than the predictor.

What a run pays is mostly the time step. In the first slab above, where a third of
the cells is weight-adjusted and steps at a tenth of what its wave speeds permit,
the clusters update 6.9 times as many cells per time as in the plain run, and the
run took 6.4 times as long; with the S impedances alone by 10, 3.9 times as many,
and for the smooth bump with every cell weight-adjusted 1.6 times as many.

Every cell of a rank carries the buffer where any of them is weight-adjusted, as
much memory again as its degrees of freedom: 1.4 kB for an elastic cell at order
4, 4.0 kB at order 6, in double precision. A rank without such a cell has none.

What is not done
----------------

- The device. Its batched kernels would take the buffer of a weight-adjusted
  cell in place of its field and apply the weight to a batch of such cells; a
  device build ignores ``weightadjustedratio`` for now.
- The media with memory variables. The anelastic block of the flux is the
  relaxation frequency times the strain rate, without material, so the memory
  variables need no weight; what is missing is a weight that leaves their rows
  alone, and, for the solver that keeps them apart, the kernels.
- The anisotropic and the poroelastic medium, whose flux a face does not carry as
  scalars per point.
- A finer choice of the cells. A cell whose material varies smoothly but by more
  than the ratio pays the time step of a jump.
- The time step at orders 2, 7 and 8, which take the smallest :math:`c` of the
  orders measured.

Checks
------

Unit checks: the operator of a medium is its inverse weight times the one without
material, in every direction, for solids, a fluid in an elastic build and
acoustic media, exactly, in Python on the declarations and in C++ against the
matrices the materials write themselves; the flux without material decomposes
into the scalars of a face as the one with the material does, between solids,
fluids and both, and at a free surface; the kernel that applies the weight does
what its matrices say where the material differs from sample to sample, a fluid
among the samples; with one material for the cell, what a weight-adjusted cell
adds -- the volume term, the local and the neighbour fluxes of a regular, a
free-surface and a fault face, the lift of a fault face -- is what a cell adds
that lifts with the material, for a solid and for a fluid; and the spread of the
impedances and the time step it gives at every order.

End to end: the runs above -- a weight-adjusted run on a constant material
against the plain one, on straight and on curved cells, elastic and acoustic; the
slabs against the plain lift and against longer time steps, at orders 3 to 6; a
smooth material against cells half as large; one rank against two.
