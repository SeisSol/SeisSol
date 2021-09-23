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

To ensure the stability of SeisSol using attenuation and local time-stepping (LTS),
it seems necessary to limit the maximum timestep using parameter ``FixTimeStep`` of the ``&Discretization`` namelist, in the main parameter file.
Else, very large elements, rarely updated by the LTS, get unstable.
Practically, a value of 0.0025 seems low enough to ensure stability.
In most cases, this should not affect the LTS speed-up.

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

   FreqCentral=2.5
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
       Qs:     return 0.1 * sqrt(mu/rho);
       Qp:     return 0.2 * sqrt(mu/rho);


FreqCentral and FreqRatio
-------------------------

| The relaxation frequencies are logarithmically equispaced, i.e.

| :math:`log(w_{i+1})-log(w_i) =` constant. 

In the parameter file, one has to give a frequency ratio of maximum to minimum frequency and a central frequency. 
For example, in the case of 3 mechanisms the following relations define the relaxation frequencies:

| :math:`w_2 = FreqCentral`  

| :math:`log(w_3)-log(w_2) = log(w_2) - log(w_1)`  

| :math:`w_3 / w_1 = FreqRatio`  

Outside of the frequency band :math:`w_1 - w_3`, Q goes to infinity, yielding
elastic behavior.

