Physical Models
===============

Overview
--------

SeisSol includes various physical models to simulate realistic earthquake scenarios.

Elastic
^^^^^^^

This is the standard model in SeisSol and it implements isotropic elastic materials. 
The constitutive behaviour is :math:`\sigma_{ij} =  \lambda \delta_{ij} \epsilon_{kk} + 2\mu \epsilon_{ij}` with stress :math:`\sigma` and strain :math:`\epsilon`. 
Elastic materials can be extended to elastoplastic materials (see :ref:`tpv-13`).

.. _anisotropic:

Anisotropic
^^^^^^^^^^^

This is an extension of the elastic material, where direction-dependent effects also play a role. 
The stress strain relation is given by the :math:`\sigma_{ij} = c_{ijkl} \epsilon_{kl}`.
Whereas isotropic materials are described by two material parameters, the tensor :math:`c` has :math:`81` entries. 
Due to symmetry considerations there are only :math:`21` independent parameters, which have to be specified in the material file:
c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66.
For more details about the anisotropic stiffness tensor, see: https://en.wikipedia.org/wiki/Hooke%27s_law#Anisotropic_materials.
All parameters have to be set, even if they are zero.
If only the two Lamé parameters are provided, SeisSol assumes isotropic behaviour. 

Anisotropy together with plasticity and dynamic rupture is not tested yet. 
You can define a dynamic rupture fault embedded in an isotropic material and have anisotropic regions elsewhere in the domain.

Viscoelastic
^^^^^^^^^^^^

Viscoelasticity is used to model the dissipation of wave energy over time. 
A full documentation can be found in :ref:`attenuation`.
