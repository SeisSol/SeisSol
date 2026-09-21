..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _damage:

Damaged material
================

What it is for
--------------

An elastic material carries a wave and is unchanged by it. A damaged material
is changed: the same strain that carries the wave also breaks the rock it
travels through, the broken rock is softer, and the softer rock carries the
next wave differently. The continuum damage-breakage rheology describes that
feedback with two internal variables and a stress that depends on them.

It is meant for what happens *around* a fault rather than on it -- the
off-fault damage zone, the softening of the medium during rupture, and the
question of how much of an earthquake's energy goes into breaking rock instead
of radiating away. On a fault itself, the friction laws of :doc:`dynamic
rupture <dynamic-rupture>` describe the same physics in a different way, and
the two are not yet joined here (see `What is supported`_).

Quantities
----------

The model is written in strain rather than stress, and it adds two scalars, so
a cell carries eleven quantities:

+-------------------------+---------------------+------------------------------------------------+
| Quantity                | Count               | Meaning                                        |
+=========================+=====================+================================================+
| :math:`\epsilon_{ij}`   | 6                   | strain, in Voigt order                         |
+-------------------------+---------------------+------------------------------------------------+
| :math:`v_i`             | 3                   | particle velocity                              |
+-------------------------+---------------------+------------------------------------------------+
| :math:`\alpha`          | 1                   | damage, between 0 and 1                        |
+-------------------------+---------------------+------------------------------------------------+
| :math:`B`               | 1                   | breakage, between 0 and 1                      |
+-------------------------+---------------------+------------------------------------------------+

Strain rather than stress is not a presentational choice. The stress of this
material is a nonlinear function of the state, so it cannot *be* the state: it
has to be evaluated, and the state has to be what it is evaluated from.
:math:`\alpha` measures how much the solid has been damaged, and
:math:`B` how much of it has lost its solid character altogether and
behaves granularly.

The equations
-------------

Write :math:`I_1 = \mathrm{tr}\,\epsilon`, :math:`I_2 = \epsilon : \epsilon`
and :math:`\xi = I_1 / \sqrt{I_2}`. The last is the *strain invariant ratio*
and it is what the model measures the kind of loading by: it is negative in
compression, positive in extension, and bounded by
:math:`\pm\sqrt{3}`.

The stress is a blend of a solid and a granular branch, weighted by the
breakage:

.. math::

   \sigma = (1 - B)\left[
       \left(\lambda_0 I_1 - \gamma_R \alpha \sqrt{I_2}\right)\mathbf{1}
       + \left(2\mu_0 - 2\gamma_R \xi_0 \alpha - \gamma_R \alpha \xi\right)\epsilon
     \right]
     + B\left[
       \left(a_1 \sqrt{I_2} + 2a_2 I_1 + 3a_3 \xi I_1\right)\mathbf{1}
       + \left(2a_0 + a_1 \xi - a_3 \xi^3\right)\epsilon
     \right]

At :math:`\alpha = 0` and :math:`B = 0` the first branch is Hooke's law with
the Lamé parameters :math:`\lambda_0` and :math:`\mu_0`, which is the check
the implementation is verified against. As :math:`\alpha` grows, the shear
modulus degrades through :math:`\gamma_R`, and the coupling term
:math:`\gamma_R \alpha \sqrt{I_2}` makes the stress depend on the *kind* of
loading and not only on its size.

The internal variables follow rate equations. Damage grows where the loading
exceeds the onset ratio :math:`\xi_0` and heals where it does not, at rates
:math:`C_d` and the healing rate; breakage grows out of damage once
:math:`\alpha` approaches a critical value, with :math:`\beta_\alpha` setting
how sharp that transition is. Both are sources rather than fluxes: they do not
propagate, they are made where the strain is.

The damage saturates at that critical value rather than at one. One is where
the variable runs out of room; the critical damage is where the solid branch
stops describing anything, and it is the smaller of the two. Past it the
tangent of the stress need no longer be positive, and a medium there carries no
waves -- a uniform state then grows roundoff at a rate that does not depend on
the timestep, because what it has lost is well-posedness rather than accuracy.
For a mild :math:`\gamma_R` the critical damage is one and the two bounds
coincide.

Momentum and the strain rate are as usual --
:math:`\rho \partial_t v = \nabla \cdot \sigma` and
:math:`\partial_t \epsilon = \tfrac{1}{2}(\nabla v + \nabla v^{T})` -- so
what makes the system nonlinear is entirely the constitutive relation and the
two sources.

Material parameters
-------------------

Per cell, from the :ref:`easi <easi>` material file:

+------------------------------------------+---------------------+------------------------+-------------------------+
| Parameter                                | SeisSol name        | Abbreviation           | Unit                    |
+==========================================+=====================+========================+=========================+
| Density                                  | ``rho``             | :math:`\rho`           | :math:`kg \cdot m^{-3}` |
+------------------------------------------+---------------------+------------------------+-------------------------+
| :math:`1^{st}` Lamé parameter, undamaged | ``lambda0``         | :math:`\lambda_0`      | :math:`Pa`              |
+------------------------------------------+---------------------+------------------------+-------------------------+
| :math:`2^{nd}` Lamé parameter, undamaged | ``mu0``             | :math:`\mu_0`          | :math:`Pa`              |
+------------------------------------------+---------------------+------------------------+-------------------------+
| Damage modulus                           | ``gammaR``          | :math:`\gamma_R`       | :math:`Pa`              |
+------------------------------------------+---------------------+------------------------+-------------------------+
| Onset strain invariant ratio             | ``xi0``             | :math:`\xi_0`          |                         |
+------------------------------------------+---------------------+------------------------+-------------------------+
| Damage rate                              | ``Cd``              | :math:`C_d`            | :math:`Pa^{-1}s^{-1}`   |
+------------------------------------------+---------------------+------------------------+-------------------------+
| Initial strain                           | ``eps_xx0`` ...     | :math:`\epsilon^0_{ij}`|                         |
|                                          | ``eps_xz0``         |                        |                         |
+------------------------------------------+---------------------+------------------------+-------------------------+

The initial strain is the background state the constitutive relation is
evaluated at: the cell's strain is added to it before the stress is formed, so
the simulated strain is a perturbation of a pre-stressed medium. The six
components are named in Voigt order, as elsewhere.

Per model, from the ``[equations]`` section of the :ref:`parameter file
<parameter-file>`:

+------------------------------------------+---------------------+------------------------+------------------+
| Parameter                                | SeisSol name        | Abbreviation           | Default          |
+==========================================+=====================+========================+==================+
| Breakage rate                            | ``breakagerate``    |                        | :math:`0`        |
+------------------------------------------+---------------------+------------------------+------------------+
| Healing rate                             | ``healingrate``     |                        | :math:`0`        |
+------------------------------------------+---------------------+------------------------+------------------+
| Width of the damage-breakage transition  | ``betaalpha``       | :math:`\beta_\alpha`   | :math:`1`        |
+------------------------------------------+---------------------+------------------------+------------------+
| Granular branch coefficients             | ``ab0`` ... ``ab3`` | :math:`a_0 \dots a_3`  | :math:`0`        |
+------------------------------------------+---------------------+------------------------+------------------+

The defaults reproduce a run with damage but without breakage or healing,
which is the configuration the model is usually first compared against. Such a
run has a definite end: the damage grows to the critical damage of the cell's
strain direction and stops there, because there is no breakage to take over.

``betaalpha`` divides and must be positive; the two rates must not be
negative. A positive breakage rate with every granular coefficient at zero is
warned about: a cell that has broken all the way then has no moduli at all, so
it carries no stress and no waves, and nothing later in a run says so. All
three are checked when the parameter file is read.

How it is solved
----------------

The material is built with the ``nonlinearck`` solver, which differs from the
linear ADER-DG scheme in three ways a user may notice.

**The stress is evaluated, not transported as a state variable.** Within a
timestep the predictor samples the cell at a set of time nodes, evaluates the
constitutive relation and the two sources at each of them, and integrates. The
time rule includes both endpoints of the step, so the internal variables march
from the start of the step to its end without a gap.

**A cell hands its neighbours more than its state.** What crosses a face is
the integrated state, the integrated stress, and two scalars from which the
dissipation of the numerical flux is scaled. The stress of a cell is a
question about that cell's material and its damage, so the cell answers it --
a neighbour never rebuilds it.

**The numerical flux carries a state-dependent speed.** Both a Rusanov and a
Godunov flux are available, as for any other material, and the dissipation of
either is scaled by the larger of the two sides' wave speeds. Those speeds
follow the damage: a softened cell carries slower waves and gets a smaller
correction. The speed is accumulated over the step rather than taken at an
instant, so it is well defined for a neighbour reading a part of the step.

The Godunov flux dissipates with the absolute value of the flux Jacobian of the
face normal. For a system whose state is a strain that is

.. math::

   |A| = \mathrm{diag}\!\left(S\,\Gamma^{-1/2}\,T,\ \Gamma^{1/2}\right),

with :math:`\Gamma` the acoustic tensor of the face normal, :math:`T` the map
from a strain to the traction it carries over the density, and :math:`S` the
lift of a vector back onto a strain. Written that way the strain rows read the
*transported stress* rather than the cell's own strain, so no material crosses
a face here either, and the form is positive semi-definite for the energy of
the system whatever pair of speeds it is given -- which is what keeps it sound
once the damage has made the tangent anisotropic.

The timestep is bounded by the stiffest state the cell can reach, over both
branches of the constitutive law and over :math:`\alpha \in [0, 1]`. A run does
not have to be restarted when a cell softens, and a granular branch stiffer
than the solid one is accounted for rather than run past.

What is supported
-----------------

+-------------------------------------------+--------------------------------------------------+
| Boundary conditions                       | regular, periodic, outflow, free surface         |
+-------------------------------------------+--------------------------------------------------+
| Local time stepping                       | yes, and as accurate as for a linear material    |
+-------------------------------------------+--------------------------------------------------+
| GPU                                       | interior faces; a face without a neighbour is    |
|                                           | refused at setup                                 |
+-------------------------------------------+--------------------------------------------------+
| Dynamic rupture                           | yes, through a nodal matrix impedance            |
+-------------------------------------------+--------------------------------------------------+
| Plasticity                                | no -- the rheology already carries the           |
|                                           | inelastic part                                   |
+-------------------------------------------+--------------------------------------------------+
| Energy output                             | momenta, kinetic energy, mean damage and         |
|                                           | breakage                                         |
+-------------------------------------------+--------------------------------------------------+
| Receivers, wavefield and surface output   | as for any other material                        |
+-------------------------------------------+--------------------------------------------------+

Three of those entries deserve their reason.

*Dynamic rupture* needs a mechanical traction at the fault, and here the
traction is derived from the state rather than being part of it: the
constitutive relation is evaluated on the fault's quadrature points and the
friction laws read the result. The impedance a node scales its Riemann problem
with is a matrix and not four scalars, because the tangent of this stress is not
isotropic -- the normal mode couples to the shear modes and the two shear modes
to each other -- and it belongs to a point and an instant rather than to the
face. In the elastic limit the rupture reproduces the linear elastic one on the
same mesh and with the same clustering to better than 0.2 percent in rupture
time, peak slip rate and final slip.

A scenario that states the fault background twice will count it twice. The
friction law compares the fault file's initial stress plus the interpolated
traction against the strength, and an initial strain in the material already
contributes to the second of those. Either the fault file or the material may
carry the background, not both.

*Local time stepping* works by reconstructing a part of a neighbour's step.
For the state that is the usual Taylor sum; for the stress it is a second
expansion, projected onto a Legendre basis in time and stored alongside. What
that reconstruction costs is the same as for a linear material: on a graded
periodic box the difference between a run with global time stepping and one at
a cluster ratio of two agrees with the linear solver's own difference to within
ten percent, and both fall at first order in the timestep. The second expansion
is therefore not what limits it -- the clustered scheme is, for either
material.

*The stored energy* is absent from the energy output rather than approximated.
The free energy of this rheology carries a square root of the second strain
invariant and terms of third order in the state; the moments the energy output
is built on reach second order in polynomials, so it cannot be assembled from
them at any accuracy. Reporting the kinetic energy and the two means, which
are exact, is preferred to a column that is named like an energy and is not
one.

Numerical notes
---------------

At :math:`I_2 = 0` -- an undeformed cell -- the ratio :math:`\xi` is
:math:`0/0`. The implementation puts a floor under :math:`I_2` that is tied to
the working precision, and sets :math:`\xi` to zero below it. An undamaged,
undeformed cell therefore behaves exactly elastically rather than producing a
NaN, which matters because it is the initial state of most runs.

The strain invariant ratio is bounded by :math:`\pm\sqrt{3}` for any strain,
and the wave-speed bound the timestep uses is taken over that range, over
:math:`\alpha \in [0, 1]`, and over both branches of the constitutive law. If a
run reports a strain invariant ratio outside those bounds, the state has left
the model's range of validity, not merely its range of accuracy.

The bound is a bound on the tangent, and a tangent is not always positive. The
critical damage the growth saturates at caps against the damage at which the
effective shear modulus vanishes, which keeps that modulus positive; the
compressional family can reach zero earlier for a strain direction that cap does
not look at. At :math:`\gamma_R = 6\times 10^{10}` and :math:`\xi = 0` it does so
at :math:`\alpha \approx 0.46`, above the critical damage that binds there. A
uniform state that will not stay uniform, at a rate independent of the timestep,
is what leaving that range looks like from the outside.

In single precision the model is more delicate than a linear one: it takes a
square root and a difference of invariants per node and timestep. Runs that
look qualitatively different between precisions should be read as a warning
about the state, not about the implementation.
