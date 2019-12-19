Physical Models
===============

Overview
--------

SeisSol includes various physical physical models to simulate realistic earthquake scenarios.

Elastic
^^^^^^^

This is the standard model in SeisSol and it implements isotropic elastic materials. 
The constitutive behaviour is :math:`\sigma =  \lambda tr(\epsilon) I + 2\mu \epsilon` with 
stress :math:`\sigma` and strain :math:`\epsilon`. Elastic materials can be extended to
elastoplastic materials (see :ref:`tpv-13`).

Anisotropic
^^^^^^^^^^^

This is an extension of the elastic material, where direction-dependent effects
also play a role. Whereas isotropic materials are described by two material parameters, the most general formulation
needs :math:`21` material parameters: :math:`\sigma_{ij} = c_{ijkl} \epsilon_{kl}`.
In the material file you can give all :math:`21` parameters giving rise to triclinic 
anisotropic material behaviour. If only the two Lamé parameters are provided, SeisSol assumes 
isotropic behaviour. Anisotropy together with plasticity and dynamic rupture is not tested yet. 
You can define a dynamic rupture fault embedded in an isotropic material and have anisotropic 
regions elsewhere in the domain.

Viscoelastic
^^^^^^^^^^^^

Viscoelasticity is used to model the dissipation of wave energy over time. 
A full documentation can be found in :ref:`attenuation`.
