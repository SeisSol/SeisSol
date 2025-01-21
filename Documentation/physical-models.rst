..
  SPDX-FileCopyrightText: 2020-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Physical Models
===============

Overview
--------

SeisSol includes various physical models to simulate realistic earthquake scenarios.
SeisSol assumes constant material values per element.
By default, the elastic and viscoelastic materials are sampled by averaging over the whole element (c.f. https://mediatum.ub.tum.de/node?id=1664043).
You can turn this feature off, by setting :code:`UseCellHomogenizedMaterial = 0` in the :ref:`parameter-file` to fall back to sampling at the element barycenter.
Anisotropic and poroelastic materials are never averaged and always sampled at the element barycenter.

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
:code:`c11`, :code:`c12`, :code:`c13`, :code:`c14`, :code:`c15`, :code:`c16`, :code:`c22`, :code:`c23`, :code:`c24`, :code:`c25`, :code:`c26`, :code:`c33`, :code:`c34`, :code:`c35`, :code:`c36`, :code:`c44`, :code:`c45`, :code:`c46`, :code:`c55`, :code:`c56`, :code:`c66`.
For more details about the anisotropic stiffness tensor, see: https://en.wikipedia.org/wiki/Hooke%27s_law#Anisotropic_materials.
All parameters have to be set, even if they are zero.
If only the two Lamé parameters are provided, SeisSol assumes isotropic behaviour.

Example of an material file for a homogeneous anisotropic fullspace. The material describes the tilted transversally isotropic medium from chapter 5.4 in
https://link.springer.com/chapter/10.1007/978-3-030-50420-5_3.

.. code:: yaml

  !ConstantMap
    map:
      rho:  2590
      c11:  6.66000000000000e10
      c12:  2.46250000000000e10
      c13:  3.44750000000000e10
      c14: -8.53035022727672e9
      c15:  0.0
      c16:  0.0
      c22:  6.29062500000000e10
      c23:  3.64187500000000e10
      c24:  4.05949408023955e9
      c25:  0.0
      c26:  0.0
      c33:  4.95562500000000e10
      c34:  7.50194506028270e9
      c35:  0.0
      c36:  0.0
      c44:  7.91875000000000e9
      c45:  0.0
      c46:  0.0
      c55:  1.40375000000000e10
      c56:  5.43430940874735e9
      c66:  2.03125000000000e10

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
