..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _attenuation:

Attenuation
===========

Introduction
------------

-  Q is measuring the decay of wave amplitude during its propagation in
   the medium and is physically defined as the ratio of wave energy to the energy dissipated per cycle of oscillation.
-  The actual dissipation mechanisms, as scattering, heat generation,
   frictional losses in the vibrating crystal lattice, etc., "causing Q"
   can be manifold and do not have to follow similar physics.
-  Q is not very frequency-dependent over a large range of (lower)
   frequencies but observations are ambiguous at higher frequencies.
-  Common choices are :math:`Q_s \sim 0.5 Q_p` and :math:`Q_s \sim 40-50V_s` (in km/s) (e.g.
   Olsen et al., 2009). :math:`Q_s = 100-200` for shallow sediments.
-  In damaged fault zones :math:`Q_s` would be expected to be as low as on the
   order of 10 and potentially smear out slip rate on the fault
   without affecting rupture speeds.

Implementation
--------------

-  Implementing attenuation in SeisSol follows the ideas of a general
   Maxwell or Zener body.
-  Each damping mechanism can be parametrized by its own relaxation
   frequency.
-  It aims at resolving a frequency-independent Q with an adequate
   number of anelastic functions and relaxation frequencies to cover the frequency range under interest.
   Usually, 3 Maxwell bodies are enough for 5% error.

Stability with Local time stepping
----------------------------------

The maximum timestep may need to be decreased in SeisSol to avoid stability issues when using attenuation and local time stepping.
Practically, we found that if the maximum timestep is below :math:`0.25 T_3` with :math:`T_3 = 1/ f_3 = 1/(\mathrm{FreqCentral} \sqrt{ \mathrm{FreqRatio}})`, stability should be ensured.
This is done by default by SeisSol when attenuation is turned on, and the parameter ``FixTimeStep`` of the the ``&Discretization`` namelist (main parameter file) is not set.
If SeisSol is yet unstable, further decrease of the maximum timestep can be tried by manually setting the value of ``FixTimeStep``.

Compiling
---------


In ccmake, use:

.. code::

    EQUATIONS                        viscoelastic2
    NUMBER_OF_MECHANISMS             3

Note that the equations='viscoelastic' is operational but deprecated.

Dispersion
----------

The attenuation implementation implies dispersion in P and S wave
velocities: that why we need to define a central frequency at which
:math:`V_p/V_s` are exact. In addition, the effective Q values are not exactly
equal to the desired Q, but are oscillating around those values. The
variation of :math:`V_p`, :math:`V_s`, :math:`Q_p` and :math:`Q_s` with frequency can be visualized using
`ViscoelasticModComp.m <https://github.com/SeisSol/SeisSol/blob/master/preprocessing/science/ViscoelasticModComp.m>`__.

Parametrisation
---------------

Add Inside the parameter file of SeisSol, in the '&equations' section
(frequencies values to be adapted to the source frequency content):

.. code:: fortran

   FreqCentral=0.5
   FreqRatio=100

The spatial variation of :math:`Q_s` and :math:`Q_p` are defined with easi in the
MaterialFileName. Here is an example of easi file, in which :math:`Q_s` and :math:`Q_p`
are directly related to the shear wave speed :math:`V_s`:

.. code:: yaml

   !ASAGI
   file: ../material/vmodel_500.nc
   parameters: [rho, mu, lambda]
   var: data
   components: !FunctionMap
     map:
       rho:    return rho;
       mu:     return mu;
       lambda: return lambda;
       Qs:     return 0.05 * sqrt(mu/rho);
       Qp:     return 0.1 * sqrt(mu/rho);


FreqCentral and FreqRatio
-------------------------

| The relaxation frequencies are logarithmically equispaced, i.e.

| :math:`\log(f_{i+1})-\log(f_i) =` constant.

In the parameter file, one has to give a frequency ratio of maximum to minimum frequency and a central frequency.
For example, in the case of 3 mechanisms the following relations define the relaxation frequencies:

| :math:`f_2 = \mathrm{FreqCentral}`

| :math:`\log(f_3)-\log(f_2) = \log(f_2) - \log(f_1)`

| :math:`f_3 / f_1 = \mathrm{FreqRatio}`

This leads  to :math:`f_1 = \mathrm{FreqCentral} / \sqrt{\mathrm{FreqRatio}}` and :math:`f_3 = \mathrm{FreqCentral}  \sqrt{\mathrm{FreqRatio}}`.

Outside of the frequency band :math:`f_1 - f_3`, Q goes to infinity, yielding
elastic behavior.

