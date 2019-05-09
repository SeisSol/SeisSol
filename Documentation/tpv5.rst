SCEC TPV5
=========

TPV5 is the first SCEC benchmark. It has spontaneous rupture on a
**vertical strike-slip fault in a homogeneous halfspace**. There are
slightly heterogeneous initial stress conditions. The earthquake rupture
is artificially nucleated in a square zone at the center of the fault
surface. The rupture then spontaneously propagates over the rest of the
fault surface. As it propagates away from the nucleation zone, it
encounters two square patches with initial stress conditions that are
different from the rest of the fault surface.

.. figure:: ./LatexFigures/tpv5_mesh.png
   :alt: Diagram of TPV5.
   :width: 9.00000cm

   Diagram of TPV5. The central square patch is the nucleation zone,
   while pink and green patches with higher and lower initial stress
   than neighbour region, respectively. 

Geometry
--------

The fault within the three-dimensional medium is a vertical
right-lateral strike-slip planar fault that resides in a half-space. The
fault reaches the Earth’s surface. The rupture is allowed within a
rectangular area that is 30000 m long :math:`\times` 15000 m deep. The
bottom boundary of and the right and left ends of the allowed 30000 m
:math:`\times` 15000 m rupture area are defined by a strength barrier.
The nucleation point is centered both along-dip and along-strike of the
30000m :math:`\times` 15000m rupture area, on the fault plane, at 15000m
along-strike and 7500m depth.

The mesh is generated in GMSH. All the files that are needed for the
simulation are provided. 

1. The tpv5.geo file contains the geometry for
the fault in a cubit region.

2. Then the .geo file can be meshed by using:

``$ gmsh tpv5.geo -3 -optimize -o tpv5.msh``

3. Then convert the .msh file to 3D Gambit neutral file

``$ gmsh2gambit -i tpv5.msh -o tpv5.neu``

The toolbox of **gmsh2gambit** is used for converting gmsh file to Gambit neutrual file. It can be found in SeisSol GitHub https://github.com/SeisSol/SeisSol/tree/master/preprocessing/meshing

4. The 3D Gambit file can be convert to PUML format for LTS in latest version of SeisSol by:
  
``$ pumgen tpv5.neu tpv5``

The compilation and usage of PUMGen can be found in https://github.com/SeisSol/PUMGen/wiki and https://seissol.readthedocs.io/en/latest/
The geometry file (.geo) can be found in this repository. 

We strongly recommend users to repeat the geometry and mesh generationg processing. However, the generated mesh file (.h5) is also available through the link (https://syncandshare.lrz.de/dl/fiNdYwqvK8cdfM5h8uRZMv9e).

.. figure:: LatexFigures/mesh5.png
   :alt: Diagram of fault geometry of TPV5. 
   :width: 10.00000cm

   Diagram of fault geometry of TPV5. The fault is 30000 m long and
   15000 m wide. The square patch has a side-length of 3000m. 

Parameters
----------

Nucleation
^^^^^^^^^^

Nucleation occurs because the initial shear stress in a 3000 m :math:`\times` 3000
m square nucleation patch is set to be higher than the initial static
yield stress in that patch. Failure occurs everywhere on the fault
plane, including in the nucleation patch, following a linear
slip-weakening fracture criterion.

TPV5 uses a linear-slip weakening friction everywhere on the fault.
There are ten parameters associated with the friction constitutive law
and fault properties in the **parameters.par**. It can be found at https://github.com/daisy20170101/SeisSol_Cookbook/.

.. literalinclude:: tpv5/parameters.par
   :language: fortran

Four friction constitutive parameters are: mu\_s, mu\_d, d\_c and
cohesion. Six stress parameters are: s\_xx, s\_yy, s\_zz, s\_xy, s\_xz,
and s\_yz. All the parameters are homogeneous on the fault except for
the nucleation patch in the center of the fault, where s\_xy is larger
compared with that elsewhere. The parameters in TPV5 are listed in Table
[table:tpv5].

+----------------------------+--------------------------------+---------+-----------------+
| Parameter                  | Description                    | Value   | Unit            |
+============================+================================+=========+=================+
| mu\_s                      | static friction coefficient    | 0.677   | dimensionless   |
+----------------------------+--------------------------------+---------+-----------------+
| mu\_d                      | dynamic friction coefficient   | 0.525   | dimensionless   |
+----------------------------+--------------------------------+---------+-----------------+
| d\_c                       | critical distance              | 0.40    | m               |
+----------------------------+--------------------------------+---------+-----------------+
| cohesion                   | friction cohesion              | 0.0     | MPa             |
+----------------------------+--------------------------------+---------+-----------------+
| s\_yy                      | stress                         | 120     | MPa             |
+----------------------------+--------------------------------+---------+-----------------+
| s\_xx,s\_zz,s\_yz,s\_xz    | stress                         | 0       | MPa             |
+----------------------------+--------------------------------+---------+-----------------+
| s\_xy                      | outside the nucleation zone    | 70      | MPa             |
+----------------------------+--------------------------------+---------+-----------------+
|                            | inside the nucleation zone     | 81.6    | MPa             |
+----------------------------+--------------------------------+---------+-----------------+

Table: Table of LSR parameters on the fault in tpv5.

Notice that there are two patches with different initial stress: the one centered at (+7.5, -7.5) has 62 MPa and (-7.5, -7.5) has 78 MPa. This inital stress is included in the fault.yaml file.

Results
~~~~~~~

All examples here can be illustrated in Paraview (Detailed instruction
can be found at ). The *output* folder contains a series of files for
fault dynamic rupture (netcdf), wave filed (netcdf), receiver (.dat) and
off-fault receivers (.dat). The fault dynamic rupture and wave filed
files can be loaded in Paraview directly. For example, open Paraview and
then go through File :math:`>>` import :math:`>>`\ prefix-fault.xdmf.

.. figure:: LatexFigures/tpv5_SRs_3s.png
   :alt: Fault slip rate in the along-strike direction
   :width: 12.00000cm

   Fault slip rate in the along-strike direction (SRs) at 4 seconds in
   TPV5, illustrated in Paraview. 

In the wave filed output file (prefix.xdmf, prefix\_vertex.h5 and
prefix\_cell.hf), the variables are shown in Table [table:wavefield]

+---------+-------------+---------------------------------+
| Index   | Parameter   | Description                     |
+=========+=============+=================================+
| 1       | U           | displacement in x-axis          |
+---------+-------------+---------------------------------+
| 2       | V           | displacement in y-axis          |
+---------+-------------+---------------------------------+
| 3       | W           | displacement in z-axis          |
+---------+-------------+---------------------------------+
| 4       | u           | particular velocity in x-axis   |
+---------+-------------+---------------------------------+
| 5       | v           | particular velocity in y-axis   |
+---------+-------------+---------------------------------+
| 6       | w           | particular velocity in z-axis   |
+---------+-------------+---------------------------------+

Table: Table of wave field output in SeisSol. Index denotes the position
used in *iOutputMask* in SeisSol parameter file.

In the fault dynamics output file (prefix-fault.xdmf,
prefix-fault\_vertex,h5 and prefix-fault\_cell,h5), the variables are
shown in Table [table:faultout]

+---------+--------------------+-------------------------------------------------------------------------------+
| Index   | Parameter          | Description                                                                   |
+=========+====================+===============================================================================+
| 1       | SRs and SRd        | slip rates in strike and dip direction                                        |
+---------+--------------------+-------------------------------------------------------------------------------+
| 2       | T\_s, T\_d, P\_n   | transient shear stress in strike and dip direction, transient normal stress   |
+---------+--------------------+-------------------------------------------------------------------------------+
| 3       | U\_n               | normal velocity (note that there is no fault opening in SeisSol)              |
+---------+--------------------+-------------------------------------------------------------------------------+
| 4       | Mud, StV           | current friction and state variable in case of RS friction                    |
+---------+--------------------+-------------------------------------------------------------------------------+
| 5       | Ts0,Td0,Pn0        | total stress, including initial stress                                        |
+---------+--------------------+-------------------------------------------------------------------------------+
| 6       | Sls and Sld        | slip in strike and dip direction                                              |
+---------+--------------------+-------------------------------------------------------------------------------+
| 7       | Vr                 | rupture velocity, computed from the spatial derivatives of the rupture time   |
+---------+--------------------+-------------------------------------------------------------------------------+
| 8       | ASl                | absolute slip                                                                 |
+---------+--------------------+-------------------------------------------------------------------------------+
| 9       | PSR                | peak slip rate                                                                |
+---------+--------------------+-------------------------------------------------------------------------------+
| 10      | RT                 | rupture time                                                                  |
+---------+--------------------+-------------------------------------------------------------------------------+
| 11      | DS                 | only with LSW, time at which ASl :math:`>` d\_c                               |
+---------+--------------------+-------------------------------------------------------------------------------+

Table: Table of fault dynamic output in SeisSol. Index denotes the
position used in *iOutputMask* in SeisSol parameter file.
