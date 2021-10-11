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

Poroelastic
^^^^^^^^^^^

In poroelastic materials a fluid and a solid phase interact with each other.
The material model introduces the pressure :math:`p` and the relative fluid velocities :math:`u_f, v_f, w_f` to the model, so that we observe 13 quantities in total.
A poroelastic material is characterised by the following material parameters:

+------------------------------------------+-----------------------+-------------------+-------------------------+
| Parameter                                |  SeisSol name         | Abbreviation      | Unit                    |
+==========================================+=======================+===================+=========================+
| Solid Bulk modulus                       |  ``bulk_solid``       | :math:`K_S`       | :math:`Pa`              |
+------------------------------------------+-----------------------+-------------------+-------------------------+
| Solid density                            |  ``rho``              | :math:`\rho_S`    | :math:`kg \cdot m^{-3}` |
+------------------------------------------+-----------------------+-------------------+-------------------------+
| Matrix :math:`1^{st}` Lamé parameter     |  ``lambda``           | :math:`\lambda_M` | :math:`Pa`              |
+------------------------------------------+-----------------------+-------------------+-------------------------+
| Matrix :math:`2^{nd}` Lamé parameter     |  ``mu``               | :math:`\mu_M`     | :math:`Pa`              |
+------------------------------------------+-----------------------+-------------------+-------------------------+
| Matrix permeability                      |  ``permeability``     | :math:`\kappa`    | :math:`m^2`             |
+------------------------------------------+-----------------------+-------------------+-------------------------+
| Matrix porosity                          |  ``porosity``         | :math:`\phi`      |                         |
+------------------------------------------+-----------------------+-------------------+-------------------------+
| Matrix tortuosity                        |  ``tortuosity``       | :math:`T`         |                         |
+------------------------------------------+-----------------------+-------------------+-------------------------+
| Fluid bulk modulus                       |  ``bulk_fluid``       | :math:`K_F`       | :math:`Pa`              |
+------------------------------------------+-----------------------+-------------------+-------------------------+
| Fluid density                            |  ``rho_fluid``        | :math:`\rho_F`    | :math:`kg \cdot m^{-3}` |
+------------------------------------------+-----------------------+-------------------+-------------------------+
| Fluid viscosity                          |  ``viscosity``        | :math:`\nu`       | :math:`Pa \cdot s`      |
+------------------------------------------+-----------------------+-------------------+-------------------------+

The implementation of poroelasticity is tested for point sources, material interfaces and free-surfaces.
Plasticity and dynamic rupture together with poroelasticity are not tested.



Viscoelastic
^^^^^^^^^^^^

Viscoelasticity is used to model the dissipation of wave energy over time. 
A full documentation can be found in :ref:`attenuation`.
