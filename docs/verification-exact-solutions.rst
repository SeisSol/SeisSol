..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Verification of dynamic rupture against exact solutions
=======================================================

The verification suite in ``tests/verification`` runs SeisSol on small meshes
that it generates itself (``minimesh.py``) and holds what the fault receivers
and the volume receivers record against solutions that are exact up to the
error of an independent integration, which is far below that of the case.
This page describes the cases that do so for dynamic rupture, the solutions
they are compared with, and how their thresholds come about.

Running the cases
-----------------

The cases are part of the tests of a build with ``TESTING=ON``. Each case is a
pair of tests per configuration, the run of SeisSol and the comparison, e.g.

.. code-block:: bash

   ctest -R "fault2d-antiplane<elastic-o6-f64>"

runs ``verification.fault2d-antiplane<elastic-o6-f64>`` and then
``verification.exact.fault2d-antiplane<elastic-o6-f64>``. All of them carry the
label ``convergence``; ``ctest -L convergence`` runs them, ``ctest -LE
convergence`` leaves them out. The meshes need ``h5py``, the spectral method
of the cases on ``fault2d-periodic`` needs ``scipy``; without it, these cases
are skipped, their runs as well as their comparisons. The cases run in the
configurations of the elastic and the anisotropic equations, the latter with an
anisotropic material that is isotropic in effect, and with media with tilted
axes of their own. With fused simulations, each simulation runs the case, and
each is compared. The cases on ``fault2d-periodic``, which take minutes each,
run in double precision without fused simulations only, and carry the label
``spectral`` as well; ``ctest -LE spectral`` leaves them out.

``fault1d.py self-test`` and ``antiplane.py self-test``
(``verification.fault1d.self-test`` and ``verification.antiplane.self-test``)
hold the exact solutions themselves against what is known in closed form or
computed otherwise, such as finite differences of a single mode.

The spectral solution of a case on ``fault2d-periodic`` is the same in every
configuration that runs it. Its comparisons keep it at the receivers in
``tests/verification/solutions`` of the build, 1.4 to 3.1 MB per case, under a
key of all that it is computed from: the inputs of the case, the versions of
numpy and scipy, and the scripts. The first comparison of a case computes it,
in one to three minutes, eight for the in-plane cases; the others, which run
one after another, and every repeated comparison read it. Removing the
directory makes the comparisons compute their solutions anew.

The meshes
----------

``fault1d``
   A column across the fault plane :math:`x = 0`, periodic along the fault
   and open across it. A fault that starts out uniform stays uniform: nothing
   is transferred along it, each side radiates a plane wave straight away
   from it, and every point of the fault solves the same ordinary
   differential equation. ``fault1d-bimaterial`` has another medium on either
   side.

``fault2d``
   The column drawn out along :math:`y` to 1200 m, for a forced rupture front
   that runs along it faster than the shear waves.

``fault2d-periodic``
   800 m along :math:`y` and :math:`266.7` m across the fault, periodic in all
   directions: what the fault sends out returns to it from the other side,
   which no outflow plane could take without reflecting part of what arrives
   inclined. A spontaneous rupture runs out on it from a nucleation along
   :math:`y`.

The exact solutions
-------------------

``fault1d.Problem``
   The 0D problem of a uniform fault: the traction is the load minus the
   radiation damping :math:`\eta V`, :math:`\eta = 1 / (1/Z_- + 1/Z_+)`, in an
   anisotropic medium a matrix, which couples the slip to the normal stress.
   It is integrated with steps far below the time step of the case, for linear
   slip weakening (FL=16), with a forced rupture time, its proxy of thermal
   pressurization (FL=1058) and its Prakash-Clifton regularization (FL=6), the
   rate-and-state laws FL=3, 4, 7 and 103, with thermal pressurization or
   without, and the imposed slip rates FL=33 and 34. The friction laws follow
   :doc:`dynamic-rupture` and :doc:`slip-rate-on-dr` and the publications they
   point to, not the code they are meant to check.

``fault1d.Front``
   A forced rupture front that runs faster than the shear waves outruns all it
   sends out: behind it, every point of the fault solves the problem of
   ``fault1d``, shifted to its rupture time, until the waves of where the front
   starts arrive.

``antiplane.Antiplane`` and ``antiplane.Inplane``
   The spectral boundary integral method of a fault between two halves of the
   same isotropic solid, periodic along :math:`y` and across the fault, for an
   antiplane rupture, which slips along :math:`z`, and for an in-plane rupture
   of mode II, which slips along :math:`y` (Perrin, Rice and Zheng, 1995;
   Geubelle and Rice, 1995; Lapusta et al., 2000). Each Fourier mode of the slip
   transfers stress along the fault through a convolution kernel, and the
   images of the fault across the period add the waves that return; in plane,
   the P waves arrive first, and the near field of P and S waves acts between
   their arrivals. Every point of the fault solves the friction law with the
   transferred stress, here with the laws of ``fault1d.Problem`` and the
   instantaneous healing of linear slip weakening (``lsw_healingThreshold``).
   The field off the fault follows from the propagators of the waves of each
   mode. The description of ``antiplane.py`` derives all kernels; its self-test
   holds them against finite differences.

The cases
---------

On ``fault1d``, the cases cover every friction law above, each with an exact
solution of its own, in isotropic and in anisotropic media; on ``fault2d``,
``fault2d-supershear`` holds a forced supershear front. On
``fault2d-periodic``:

=========================================== ====== ==========================================
Case                                        FL     What it adds
=========================================== ====== ==========================================
``fault2d-antiplane``                       16     a spontaneous antiplane rupture
``fault2d-antiplane-rs-*``                  3-103  the rate-and-state laws
``fault2d-antiplane-rs-fvw-tp``             103    thermal pressurization
``fault2d-antiplane-lsw-bimaterial``        6      the Prakash-Clifton regularization
``fault2d-antiplane-lsw-tpapprox``          1058   the proxy of thermal pressurization
``fault2d-antiplane-lsw-healing``           16     instantaneous healing
``fault2d-antiplane-imposed-*``             33, 34 an imposed slip that runs out along y
``fault2d-inplane``                         16     a spontaneous in-plane rupture
``fault2d-inplane-supershear``              16     its transition to supershear
=========================================== ====== ==========================================

The imposed slip rate FL=35 imposes all of the slip within the sub-step of the
onset: its slip rate, and the field it sends out, depend on the time step and
converge to no function; it has no comparison. A fault between two different
media in plane, where the slip changes the normal stress, has none either: the
spectral method takes one medium on both sides, and on a mesh periodic across
the fault, the second interface between the media, halfway to the next fault,
reflects what arrives there.

What is compared
----------------

The fault receivers report the slip rate, the normal stress and the direction
of the slip at their own position, from the velocities there; the other
quantities come from the Gauss point of the fault nearest to them: the
traction, which takes the friction coefficient and the initial stress of the
Gauss point, the slip and the accumulated slip, the peak slip rate, the
friction coefficient and the state variable, the pore pressure and the
temperature of thermal pressurization, and the rupture time and the dynamic
stress time (see :doc:`fault-output`). On ``fault2d-periodic``, the fault
receivers lie on Gauss points of the fault of the quadrature of order 6, so
that all of this holds at their position; the comparison checks that they do,
for the quadrature of the build, and compares the quantities of the Gauss
points only where they do. Where none of them does, as with Stroud's rules,
whose points depend on the order of the vertices of a face, the thresholds of
those quantities are not applied either.

Each quantity is compared as its largest deviation over time, relative to the
largest value of the exact solution, the events as their deviation relative
to their time. Where the receivers report something else than the physical
quantity, they are compared with what they report: the traction of an imposed
slip rate is that of the Godunov state, i.e. of the welded fault; with thermal
pressurization, the normal stress is the effective one, and the pore pressure
a stress, negative where it pushes the faces apart. Where slip weakening
heals, the slip rate, the traction and the friction coefficient jump; the time
it heals, where the friction coefficient first rises, is compared as an event,
and the jumps themselves are left out. The jump of the slip rate runs along the
fault faster than shear waves, and sends out a jump of the field, which the
method resolves only over the time the healing takes through a cell; the volume
receivers are compared up to the time the fault first heals.

Thresholds
----------

``tests/verification/cases/thresholds-fault1d.json`` holds a threshold per
case, quantity and key ``case/equations/o<order>/<precision>``. A key without
thresholds passes, and prints the deviations. The thresholds are 1.5 times the
largest deviation over the configurations of the key, rounded up to two digits,
and no less than a millionth for the direction and the peak slip rate, which
may deviate by no more than the rounding: they catch a change of the numerics,
not only a gross error. On ``fault2d-periodic``, the four builds of order 6,
elastic and anisotropic, in double and in single precision, deviated alike when
the thresholds were taken, to a few per cent where the deviations exceed the
rounding; the two builds in double precision, which run these cases, share the
thresholds of the largest deviation among the four. The comparison prints its
deviations in the form of the table with ``--record``, or with
``SEISSOL_RECORD_THRESHOLDS`` set; a change that moves them has to say so in
the thresholds.
