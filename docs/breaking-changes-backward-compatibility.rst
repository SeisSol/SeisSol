..
  SPDX-FileCopyrightText: 2022 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Breaking changes in backward compatibility
==========================================

To keep up-to-date with changes in compute-centers and geoscientists' needs, breaking changes sometimes needed to be introduced in SeisSol.

All breaking changes for version 0.9.0 and later are listed here.

The fault tag of the elementwise fault output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(since the unification of the output modules)

The ``fault-tag`` cell field of the elementwise fault output carried the
identifier of the face rather than the tag the mesh gave it, which is what the
``global-id`` field beside it holds, so the files had the identifier twice and
the tag not at all. Post-processing that read ``fault-tag`` and got what it
expected was reading an identifier; one that grouped by it was grouping by face.

The on-fault and off-fault receiver files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(since the unification of the output modules)

``receiverFormat = 'hdf5'`` writes a different file. It used to be one wide
table with a column count taken from the widest receiver, described by the
attributes ``DimNames`` and ``VariableNames``; it is now a dataset per quantity
set, with the quantities as the members of a compound and the receivers
described by columns beside it. :ref:`off_fault_receivers` has the layout.
Post-processing of ``-receivers.h5`` has to follow.

The text output, which is the default, is unchanged, and so is the file name and
the write interval. The on-fault receivers gained the same HDF5 layout as an
option under ``format = 'hdf5'`` in the ``Pickpoint`` section, but keep writing
text unless it is asked for.

The tables beside the output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(since the unification of the output modules)

The column names of ``-clustering.csv``, ``-threadPinning.csv`` and
``-miniSeissol.csv`` are now quoted, like every other piece of text in those
files. What the columns are and what they hold is unchanged. A reader using a
CSV parser needs no change; one comparing the header line literally does.

Output file names
~~~~~~~~~~~~~~~~~
(since the unification of the output modules)

Every mesh output is now named after what it holds, so three file names changed.
Post-processing that opens them by name has to follow.

* The wavefield written through Xdmf was ``<prefix>.xdmf``; it is now
  ``<prefix>-wavefield.xdmf``, next to the ``<prefix>-wavefield.vtkhdf`` that the
  high-order output already used. The bare prefix carried no indication of what
  was in the file.
* The high-order free-surface output was ``<prefix>-free-surface.vtkhdf``; it is
  now ``<prefix>-surface.vtkhdf``, which is what the Xdmf free-surface output has
  always been called.
* The high-order elementwise fault output was ``<prefix>-fault-elementwise.vtkhdf``;
  it is now ``<prefix>-fault.vtkhdf``, matching the Xdmf fault output. The output
  it is distinguished from -- the on-fault receivers -- is written by a different
  module under a different name, so the qualifier distinguished nothing.

Refined wavefield output
~~~~~~~~~~~~~~~~~~~~~~~~
(since the unification of the output modules)

Two things changed about the volume output at the same time, and both affect a comparison against
files written by an older version.

First, a refined wavefield output used to sample the subcells in a different vertex labelling than
the one the output mesh was built with, so the value written for a subcell was the solution at the
center of one of its siblings. With ``refinement = 1``, three of the four subcells of every element
carried a neighbour's value; with ``refinement = 2`` and ``3``, the inner subcells were sampled at
a point that is not the center of any subcell at all. ``refinement = 0`` was unaffected, because
the center of the whole element does not move under that relabelling. This is now corrected, so
output written with ``refinement > 0`` differs from what older versions produced.

Second, how the solution reaches the output points is now a parameter rather than a property of
the writer. It is ``wavefieldprojection`` for the wavefield and ``surfaceprojection`` for the free
surface (see :ref:`wave_field_output` and :ref:`free_surface_output`). The defaults reproduce what
each output did before -- ``pointwise`` for the wavefield, ``l2`` for the free surface -- so no
parameter file needs to change; the option exists to make the two comparable when that is wanted.

Energy Output
~~~~~~~~~~~~~
(since 0.9.0, `#531 <https://github.com/SeisSol/SeisSol/pull/531>`_, April 2022)

Since we merged GitHub pull request `#531 <https://github.com/SeisSol/SeisSol/pull/531>`_ (April 2022), the seismic moment time history output,
from which the moment rate can be post-processed, is integrated into the energy output (see :ref:`energy_output`).
Therefore, the parameters `magnitude_output_on`, `energy_rate_output_on` and `energy_rate_printtimeinterval` have been removed from the `DynamicRupture` namelist in the main parameter file.

C++ dynamic rupture implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(1.0.0 to 1.3.2, since `#625 <https://github.com/SeisSol/SeisSol/pull/625>`_, September 2022)

While porting dynamic rupture to C++, we changed a few parameter names to make things more consistent.
The new dynamic rupture implementation has been merged in September 2022 (GitHub pull request `#625 <https://github.com/SeisSol/SeisSol/pull/625>`_).
The linear slip weakening friction laws FL=2 (nucleation by stress increase) and FL=16 (forced time rupture nucleation) have been merged (the new friction law is FL=16; but FL=2 will be accepted as input after `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_).
Because of this change, FL=16 now requires nucleation stress or tractions to be specified in the fault-specific yaml file.

Parameter file (`parameters.par`):

+-----------------+-------------------------------------------------------------------------------------------------+
| old             | new                                                                                             |
+=================+=================================================================================================+
| ``0d0``         | ``0.0`` (only until `#1173 <https://github.com/SeisSol/SeisSol/pull/1173>`_)                    |
+-----------------+-------------------------------------------------------------------------------------------------+
| ``v_star``      | ``pc_vStar`` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)               |
+-----------------+-------------------------------------------------------------------------------------------------+
| ``L``           | ``pc_prakashLength`` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)       |
+-----------------+-------------------------------------------------------------------------------------------------+
| ``mu_w``        | ``rs_muW`` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)                 |
+-----------------+-------------------------------------------------------------------------------------------------+
| ``alpha_th``    | ``tp_thermalDiffusivity`` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)  |
+-----------------+-------------------------------------------------------------------------------------------------+
| ``rho_c``       | ``tp_heatCapacity`` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)        |
+-----------------+-------------------------------------------------------------------------------------------------+
| ``tp_lambda``   | ``tp_undrainedTPResponse`` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_) |
+-----------------+-------------------------------------------------------------------------------------------------+
| ``initemp``     | ``tp_iniTemp`` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)             |
+-----------------+-------------------------------------------------------------------------------------------------+
| ``inipressure`` | ``tp_iniPressure`` (only until `#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_)         |
+-----------------+-------------------------------------------------------------------------------------------------+

Fault-specific yaml file (`fault.yaml`):

+-------------------------------+---------------------------------------------------------------------------------------------------+
| old                           | new                                                                                               |
+===============================+===================================================================================================+
| ``RS_sl0``                    |  ``rs_sl0`` (only until `#1427 <https://github.com/SeisSol/SeisSol/pull/1427>`_)                  |
+-------------------------------+---------------------------------------------------------------------------------------------------+
| ``alpha_hy``                  |  ``tp_hydraulicDiffusivity`` (only until `#1427 <https://github.com/SeisSol/SeisSol/pull/1427>`_) |
+-------------------------------+---------------------------------------------------------------------------------------------------+
| ``TP_half_width_shear_zone``  |  ``tp_halfWidthShearZone`` (only until `#1427 <https://github.com/SeisSol/SeisSol/pull/1427>`_)   |
+-------------------------------+---------------------------------------------------------------------------------------------------+
| ``Ts0``                       |  ``T_s`` (only until `#1427 <https://github.com/SeisSol/SeisSol/pull/1427>`_)                     |
+-------------------------------+---------------------------------------------------------------------------------------------------+
| ``Td0``                       |  ``T_d`` (only until `#1427 <https://github.com/SeisSol/SeisSol/pull/1427>`_)                     |
+-------------------------------+---------------------------------------------------------------------------------------------------+
| ``Pn0``                       |  ``T_n`` (only until `#1427 <https://github.com/SeisSol/SeisSol/pull/1427>`_)                     |
+-------------------------------+---------------------------------------------------------------------------------------------------+
| ``RS_f0``                     |  ``rs_f0`` (only until `#1540 <https://github.com/SeisSol/SeisSol/pull/1540>`_)                   |
+-------------------------------+---------------------------------------------------------------------------------------------------+
| ``RS_b``                      |  ``rs_b`` (only until `#1540 <https://github.com/SeisSol/SeisSol/pull/1540>`_)                    |
+-------------------------------+---------------------------------------------------------------------------------------------------+
| ``RS_muw``                    |  ``rs_muw`` (only until `#1540 <https://github.com/SeisSol/SeisSol/pull/1540>`_)                  |
+-------------------------------+---------------------------------------------------------------------------------------------------+

FORTRAN Removal
~~~~~~~~~~~~~~~

(since 1.1.0, `#829 <https://github.com/SeisSol/SeisSol/pull/829>`_, April 2023)

All FORTRAN had been removed; and several parameters have been deprecated.
However, all previous configuration files continue working.

Relative Paths
~~~~~~~~~~~~~~

(since 1.2.0, `#1156 <https://github.com/SeisSol/SeisSol/pull/1156>`_, August 2024)

SeisSol now looks for all paths in the parameter file first relatively to the parameter file,
and only then relative to the execution directory. This rule does not apply for the output directory.

Before this change, also all other paths were only viewed relative to the execution directory.

For absolute paths, nothing has changed.

DR Traction Computation
~~~~~~~~~~~~~~~~~~~~~~~

(since 1.1.0, `#895 <https://github.com/SeisSol/SeisSol/pull/895>`_, September 2023)

There was an adjustment in the Rate-and-State friction law computation; thus the results differ slightly between 1.0.1, and 1.1.0 onwards for the CPU computation.
The adjustment was propagated to the GPU implementation after 1.3.0.

The respective `commit <https://github.com/SeisSol/SeisSol/commit/73b284b7a8a2323170766f3ab594312a31f514c1>`_
updated the CPU implementation; the GPU implementation was updated by
`#1288 <https://github.com/SeisSol/SeisSol/pull/1288>`_.

Name of the Strain Rate Output for Off-Fault Receivers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(since 1.2.0, `#1126 <https://github.com/SeisSol/SeisSol/pull/1126>`_, June 2024)

The strain rate output was named just "strain" output for the off-fault receivers.
The corresponding option was likewise called :code:`ReceiverComputeStrain`,
not :code:`ReceiverComputeStrainRate`.

Poroelastic Time Basis
~~~~~~~~~~~~~~~~~~~~~~

(unreleased, `#1374 <https://github.com/SeisSol/SeisSol/pull/1374>`_, August 2025)

Poroelasticity is solved with a space-time predictor, whose coefficients are expressed in a Legendre
basis in time, while all other materials use a monomial (Taylor) basis.

Up to this fix, two areas were evaluated
with the same (monomial) basis on all cases, affecting LTS and the fault stress output.

Poroelastic results may therefore differ from earlier versions.

Poroelastic Dynamic Rupture Impedance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(unreleased)

The shear impedance of a poroelastic fault is now :math:`Z_s = \sqrt{\mu \rho_1}`, with the
statically condensed density :math:`\rho_1 = \bar\rho - \rho_f^2 / m` that the Biot system
propagates shear waves with, instead of the density of the solid grains. The scalar impedances used
by the friction update, the slip accumulation and the fault receiver output now come from the same
matrix as the Riemann solver, and the wave impedance itself is computed in closed form rather than
from an eigendecomposition.

Results of poroelastic dynamic rupture simulations change accordingly; how much depends on the
porosity and the tortuosity. For the material values of the poroelastic test cases the impedance
drops by 5 to 15 percent.

Frictional Energy of a Bimaterial Fault
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(unreleased)

The interface traction the frictional energy is integrated with is :math:`\tau^* = b^+ \tau^+ +
b^- \tau^-`, with :math:`b^\pm = \eta Y^\pm`. The energy computation paired :math:`b^+` with the
traction of the *minus* side and vice versa, which disagreed with both the Riemann solver and the
:code:`computeTractionInterpolated` kernel the fault output uses.

The two coefficients are equal for a fault with the same material on both sides, so only
bimaterial faults are affected. The frictional energy in the energy output changes there; the
simulation itself does not.

Anisotropic Eigenbasis
~~~~~~~~~~~~~~~~~~~~~~

(unreleased)

One entry of the eigenbasis an anisotropic material is transformed with carried :math:`c_{46}`
where the derivation asks for :math:`c_{56}`, i.e. the coupling of :math:`\sigma_{xy}` to
:math:`u_z` instead of the intended one.

The matrix enters the boundary conditions, so anisotropic simulations with a free surface or an
absorbing boundary change. Materials for which both coefficients vanish -- isotropy, VTI with the
symmetry axis along a coordinate axis -- are unaffected.

Potency and Seismic Moment Quadrature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(since before 0.9.0, `#527 <https://github.com/SeisSol/SeisSol/pull/527>`_, April 2022)

The potency and the seismic moment were computed by averaging the value over all points.
Now, to make the integration more exact, they are instead weighed by the quadrature rule
the underlying DR implementation uses. As a result, the computed seismic moment and magnitude
may slightly change compared to before.
