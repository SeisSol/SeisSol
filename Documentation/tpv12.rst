SCEC TPV12
==========

TPV12 and 13 are recommended by SCEC for elastic/plastic wave
propagation code validation. TPV 12 describes spontaneous rupture on a
**60-degree dipping normal fault** in a homogeneous half-space. Material
properties are linear **elastic**. Initial stress conditions are
dependent on depth. Strongly super-shear rupture conditions.

.. figure:: LatexFigures/tpv12_13.png
   :alt: Diagram of geometry of TPV12/13.
   :width: 9.00000cm

   Diagram of geometry of TPV12/13. 60-degree dipping normal fault.

Geometry
~~~~~~~~

The model volume is a half-space. The fault is a 60-degree dipping,
planar, normal fault. The fault reaches the Earth’s surface. Rupture is
allowed within a rectangular area measuring 30000 m along-strike and
15000 m down-dip.

Note that 15000 m down-dip corresponds to a depth of 12990.38 m. A node
which lies exactly on the border of the 30000 m :math:`\times` 15000 m
rectangle is considered to be inside the rectangle, and so should be
permitted to rupture.

The portions of the fault below, to the left of, and to the right of the
30000 m :math:`\times` 15000 m rectangle are a strength barrier, within
which the fault is not allowed to rupture.

The nucleation zone is a square measuring 3000 m × 3000 m. The center of
the square is located 12000 m down-dip (at a depth of 10392.30 m), and
is centered along-strike.

The geometry is generated with GMSH. All the files that are needed for
the simulation are provided in

.. figure:: LatexFigures/tpv12mesh2.png
   :alt: Diagram of a 60-degree dipping fault in Gmsh.
   :width: 9.00000cm

   Diagram of a 60-degree dipping fault in Gmsh. The surrouding box is
   500 km long and 500 km wide and 50 km hight. The fault cuts through
   the free surface. 

The geometry and mesh generation process is similar as TPV5. The
planar-fault geometry is build with Gmsh (Figure [fig:tpv12geo]). All
the files that are needed for the simulation are provided in .

Nucleation
~~~~~~~~~~

In previous benchmarks, nucleation was achieved by imposing a higher
initial shear stress within a nucleation zone. In TPV12 and TPV13,
nucleation is achieved by selecting a lower static coefficient of
friction within a nucleation zone, so that the initial shear stress
(which is implied by the initial stress tensor) is greater than the
yield stress.

Outside the 30000 m \* 15000 m rectangular rupture area there is a
strength barrier, where nodes are not allowed to slip. Some codes
implement the strength barrier by setting the static coefficient of
friction and frictional cohesion to very large values. Other codes
implement the strength barrier in other ways.

Parameters
~~~~~~~~~~

LSR parameters
^^^^^^^^^^^^^^

TPV12 uses a linear slip weakening law on the fault with different
parameters inside and outside the nucleation zone. The parameters are
listed in Table [table:tpv12lsr].

+-------------+--------------------------------+------------+--------+
| Parameter   | inside the nucleation zone     | Value      | Unit   |
+=============+================================+============+========+
|               inside the nucleation zone                           |
+-------------+--------------------------------+------------+--------+
| mu\_s       | static friction coefficient    | 0.54       |        |
+-------------+--------------------------------+------------+--------+
| mu\_d       | dynamic friction coefficient   | 0.10       |        |
+-------------+--------------------------------+------------+--------+
| d\_c        | critical distance              | 0.50       | m      |
+-------------+--------------------------------+------------+--------+
| cohesion    | shear stress cohesion          | -200 000   | Pa     |
+-------------+--------------------------------+------------+--------+
|               outside the nucleation zone                          |
+-------------+--------------------------------+------------+--------+
| mu\_s       | static friction coefficient    | 0.70       |        |
+-------------+--------------------------------+------------+--------+
| mu\_d       | dynamic friction coefficient   | 0.10       |        |
+-------------+--------------------------------+------------+--------+
| d\_c        | critical distance              | 0.50       | m      |
+-------------+--------------------------------+------------+--------+
| cohesion    | shear stress cohesion          | -200 000   | Pa     |
+-------------+--------------------------------+------------+--------+

Table: Table of LSR parameters on the fault.

Initial stress
~~~~~~~~~~~~~~

The initial stress on the fault is depth-dependent in TPV12/13. In the
shallower portion above 11951.15 m, the stress field is optimal
orientated while the other is isotropic.

+-----------------------------------+----------------------------------+
|   Parameter                       |       Value                      |
+===================================+==================================+
|   above 11951.15 m                                                   |
+-----------------------------------+----------------------------------+
| :math:`\sigma_1`                  |  26460 Pa/m * H                  |
+-----------------------------------+----------------------------------+
| :math:`\sigma_3`                  |  15624.3 Pa/m * H                |
+-----------------------------------+----------------------------------+
| :math:`\sigma_2`                  |  :math:`(\sigma_1+\sigma_3)/2`   |
+-----------------------------------+----------------------------------+
| :math:`P_f`                       | :math:`1000 kg/m^3 *9.8 m/s^2 *H`|
+-----------------------------------+----------------------------------+
|   below 11951.15 m                                                   |
+-----------------------------------+----------------------------------+
|:math:`\sigma_1,\sigma_2,\sigma_3` | :math:`2700 kg/m^3 *9.8 m/s^2 *H`|
+-----------------------------------+----------------------------------+


Results
~~~~~~~

SeisSol output xdmf file that can be loaded in Paraview directly. The
wave field and fault output files have the same format as in TPV5.

.. figure:: LatexFigures/SR_W_tpv12.png
   :alt: Paraivew figure of TPV12 output.
   :width: 11.00000cm

   Paraivew figure of TPV12 output. Fault slip rate in dip-direction
   (SRd) and vertical velocity (w) in wave field. The roughed cutoff
   surface demonstrates the unstructured tetrahedral meshing. 
   
