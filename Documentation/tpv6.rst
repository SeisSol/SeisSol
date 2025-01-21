..
  SPDX-FileCopyrightText: 2019-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _tpv6:

SCEC TPV6
=========

TPV6 is intended to reside in the “well-posed” regime for **bimaterial problems** and so uses a very high shear modulus (density\*vs\*vs) contrast. Material properties are homogeneous within each side of the fault, but change when one traverses to the other side of the fault.

+----------------+-------------------------------+---------+---------+
| Parameter      | Description                   | Value   | Unit    |
+================+===============================+=========+=========+
| Vp             | P velocity of the far side    | 3750    | m/s     |
+----------------+-------------------------------+---------+---------+
| Vs             | S velocity of the far side    | 2165    | m/s     |
+----------------+-------------------------------+---------+---------+
| :math:`\rho`   | density of the far side       | 225     | kg/m3   |
+----------------+-------------------------------+---------+---------+
| Vp             | P velocity of the near side   | 6000    | m/s     |
+----------------+-------------------------------+---------+---------+
| Vs             | S velocity of the near side   | 3464    | m/s     |
+----------------+-------------------------------+---------+---------+
| :math:`\rho`   | density of the near side      | 2670    | kg/m3   |
+----------------+-------------------------------+---------+---------+

Table: Table of bi-material parameters in Tpv6.

Geometry
~~~~~~~~

TPV6 uses the same geometry as TPV5. The fault within the
three-dimensional medium is a vertical right-lateral strike-slip planar
fault that resides in a half-space. The fault reaches the Earth’s
surface. The rupture is allowed within a rectangular area that is 30000
m long :math:`\times` 15000 m deep. The bottom boundary of and the right
and left ends of the allowed 30000 m :math:`\times` 15000 m rupture area
are defined by a strength barrier. The nucleation point is centered both
along-dip and along-strike of the 30000m :math:`\times` 15000m rupture
area, on the fault plane, at 15000m along-strike and 7500m depth.

Parameters
~~~~~~~~~~

TPV6 uses a similar parameter setup as TPV5 except for the bulk
parameters.

+-------------------+-----------------------------------+---------+---------+
| Parameter         | Description                       | Value   | Unit    |
+===================+===================================+=========+=========+
| :math:`\rho`      | density of the far side           | 2225    | kg/m3   |
+-------------------+-----------------------------------+---------+---------+
| :math:`\lambda`   | Lame parameter of the far side    | 10.4    | GPa     |
+-------------------+-----------------------------------+---------+---------+
| :math:`\mu`       | Lame parameter of the far side    | 10.4    | GPa     |
+-------------------+-----------------------------------+---------+---------+
| :math:`\rho`      | density of the near side          | 2670    | kg/m3   |
+-------------------+-----------------------------------+---------+---------+
| :math:`\lambda`   | Lame parameter of the near side   | 32      | GPa     |
+-------------------+-----------------------------------+---------+---------+
| :math:`\mu`       | Lame parameter of the near side   | 32      | GPa     |
+-------------------+-----------------------------------+---------+---------+

Table: Table of bi-material parameters used in SeisSol for Tpv6.

Results
~~~~~~~

Figure [fig:tpv6-4s] and [fig:tpv6-7s] show the fault slip rate at 4 s
and 7 s, respectively. The slip front is asymmetric when compared with
TPV5 (Figure [fig:tpv5-4s]). Figure [fig:tpv6\_velocity] shows velocity
recorded at two off-fault receivers. The wave picks arrives at the
far-side receiver lower than those at the near-side receiver.

.. figure:: LatexFigures/tpv6_SRs_4s.jpg
   :alt: Fault slip rate at 4 seconds in the along-strike direction in
   :width: 12.00000cm
   :align: center

   Fault slip rate at 4 seconds in the along-strike direction in TPV6.

.. figure:: LatexFigures/tpv6_SRs_7s.jpg
   :alt: Fault slip rate at 7 seconds in the along-strike direction
   :width: 12.00000cm
   :align: center

   Fault slip rate at 7 seconds in the along-strike direction in TPV6.

.. figure:: LatexFigures/tpv6_off_velocity.png
   :alt: Particle velocity at two opposite stations across the fault (+/- 9 km normal to the fault).
   :width: 12.00000cm
   :align: center
