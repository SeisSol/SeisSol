..
  SPDX-FileCopyrightText: 2019-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _tpv12:

SCEC TPV12
==========

TPV12 and 13 are verification benchmarks for  a 60 degree dipping normal fault embedded in a homogeneous halfspace.
The rheology is linear elastic in TPV12 and non-associative Drucker-Prager visco-plastic in TPV13.
Initial stress conditions are
dependent on depth. Strongly super-shear rupture conditions are assumed.

.. figure:: LatexFigures/tpv12_13.png
   :alt: Geometry in SCEC benchmarks TPV12/13
   :width: 9.00000cm
   :align: center

   Figure 1: Geometry in SCEC benchmarks TPV12/13, including a 60-degree dipping normal fault.

Geometry
~~~~~~~~

A half-space domain is assumed. The fault is a 60-degree dipping, 30 km :math:`\times` 15 km
planar, normal fault. It reaches the Earth’s surface, and extends up to 12.99038 km depth.
A square nucleation area of 3 :math:`\times` 3 km is assumed.
In fault coordinates, it is centered at (0, -12) km, that is at a depth of 10.3923 km.
The CAD model (see :ref:`Figure 2 <TP12_figure_2>`) and mesh are generated with `Gmsh <https://gmsh.info/>`_. All the files that are needed for
the simulation are provided `here <https://github.com/SeisSol/Examples/tree/master/tpv12_13>`_.
The geometry and mesh generation process is similar to TPV5.

.. _TP12_figure_2:

.. figure:: LatexFigures/tpv12mesh2.png
   :alt: TPV12/13 geometry modeled in Gmsh.
   :width: 9.00000cm
   :align: center

   Figure 2: TPV12/13 geometry modeled in Gmsh. The domain box is
   500 km :math:`\times` 500 km :math:`\times` 50 km. The fault reaches the top surface.


Nucleation strategy
~~~~~~~~~~~~~~~~~~~

In previous benchmarks, nucleation is achieved by imposing higher
initial shear stress within the nucleation zone. In TPV12 and TPV13,
nucleation is achieved by decreasing the static friction coefficient,
leading the initial shear stress to exceed fault strength.


Parameters
~~~~~~~~~~

friction parameters
^^^^^^^^^^^^^^^^^^^^

TPV12 uses a linear slip weakening law with different
parameters inside and outside the nucleation zone. The parameters are
listed in the table below:

+-------------+--------------------------------+------------+--------+
| Parameter   | description                    | Value      | Unit   |
+=============+================================+============+========+
| d\_c        | critical distance              | 0.50       | m      |
+-------------+--------------------------------+------------+--------+
| cohesion    | shear stress cohesion          | -200 000   | Pa     |
+-------------+--------------------------------+------------+--------+
|               inside the nucleation zone                           |
+-------------+--------------------------------+------------+--------+
| mu\_s       | static friction coefficient    | 0.54       |        |
+-------------+--------------------------------+------------+--------+
| mu\_d       | dynamic friction coefficient   | 0.10       |        |
+-------------+--------------------------------+------------+--------+
|               outside the nucleation zone                          |
+-------------+--------------------------------+------------+--------+
| mu\_s       | static friction coefficient    | 0.70       |        |
+-------------+--------------------------------+------------+--------+
| mu\_d       | dynamic friction coefficient   | 0.10       |        |
+-------------+--------------------------------+------------+--------+

Table 1: frictions parameters used in TPV12/13.

Initial stress
~~~~~~~~~~~~~~

The initial stress is assumed depth-dependent in TPV12/13. Above 11951.15 m depth, deviatoric stresses are non-zero, and the stress field is well
orientated for rupture of the normal fault. Below, the stress is isotropic.

+-----------------------------------+-----------------------------------------------------------------------+
|   Parameter                       |       Value                                                           |
+===================================+=======================================================================+
|   above 11951.15 m depth                                                                                  |
+-----------------------------------+-----------------------------------------------------------------------+
| :math:`\sigma_1`                  |  26460 Pa/m :math:`\times` H                                          |
+-----------------------------------+-----------------------------------------------------------------------+
| :math:`\sigma_3`                  |  15624.3 Pa/m :math:`\times` H                                        |
+-----------------------------------+-----------------------------------------------------------------------+
| :math:`\sigma_2`                  |  :math:`(\sigma_1+\sigma_3)/2`                                        |
+-----------------------------------+-----------------------------------------------------------------------+
| :math:`P_f`                       | :math:`1000 \mathrm{kg/m}^3 \times 9.8 \mathrm{m/s}^2 \times H`       |
+-----------------------------------+-----------------------------------------------------------------------+
|   below 11951.15 m depth                                                                                  |
+-----------------------------------+-----------------------------------------------------------------------+
| :math:`\sigma_1,\sigma_2,\sigma_3`| :math:`2700 \mathrm{kg/m}^3 \times 9.8 \mathrm{m/s}^2 \times H`       |
+-----------------------------------+-----------------------------------------------------------------------+


Results
~~~~~~~

SeisSol output can be visualized directly in Paraview by loading their xdmf files.

.. figure:: LatexFigures/SR_W_tpv12.png
   :alt: fault and volume output of TPV12 vizualized in Paraview.
   :width: 11.00000cm
   :align: center

   Figure 3: Fault and volume output of TPV12 visualized in Paraview. Fault slip rate in dip-direction
   (SRd) and vertical velocity (w) in the volume. A cut-view of the volume output allows visualizing the unstructured tetrahedral mesh.
