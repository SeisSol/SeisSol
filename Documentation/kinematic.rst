Kinematic source example - 1994 Northridge earthquake
=====================================================

We use this earthquake to demonstrate how to setup dynamic rupture model
with kinematic rupture source in SeisSol.

The 1994 Northridge earthquake occurred on January 17, at 4:30:55 a.m.
PST and had its epicenter in Reseda, a neighborhood in the north-central
San Fernando Valley region of Los Angeles, California, USA. It had a
duration of approximately 10–20 seconds. The blind thrust earthquake had
a magnitude of 6.7 (Mw). This is a typical reverse-slip earthquake. The
fault orients to N122\ :math:`^\circ`\ E and dips at 40\ :math:`^\circ`.
The simulation can be used to build similar model with moderate
modifications.

Geometry
~~~~~~~~

The fault geometry is made in Gmsh. Fault: plane fault 20 km\*25 km
dipping at 40-degree.

Region: 100 km\*100 km \*60 km.

.. figure:: LatexFigures/1994northridge.png
   :alt: Geometry of 1994 northridge earthquake.
   :width: 12.00000cm

   Geometry of 1994 northridge earthquake. A planar fault orients at 122
   degree and dip at 40 degree. The dimension of fault is 20 km along
   strike and 25 km along down-dip.

Kinematic rupture Source
~~~~~~~~~~~~~~~~~~~~~~~~

The kinematic source of the earthquake can be found at . The *standard
rupture format* can be used directly in SeisSol, with the following
lines in parameter.par file.

::
  
  &SourceType
  Type = 42
  FileName=’northridge.nrf’
  /

Download standard rupture format file (northridge.srf) can be found in https://scec.usc.edu/scecpedia/Standard_Rupture_Format.
Please note that the SCEC units are different with SeisSol units in some
aspect.

The fault are divided in to 20 grids along the strike and 25 grids
  along the dip. The source time function (STF) of each rectangular
  elements is given in the file , whose format looks like the following:
  
::

  verison (1.0)
  PLANE 1
  ELON ELAT NSTK NDIP LEN WID STK DIP DTOP SHYP DHYP
  POINTS 500
  LON LAT DEP STK DIP AREA TINIT DT
  RAKE SLIP1 NT1 SLIP2 NT2 SLIP3 NT3
  SR1[1] SR1[2] SR1[3] . . . SR1[NT1]
  SR2[1] SR2[2] SR2[3] . . . SR2[NT3]
  SR3[1] SR3[2] SR3[3] . . . SR3[NT3]
  ... 

Explanations for the input file:

**Line 1**: version

**Line 2**: Number of fault planes

**Line 3**:
ELON top center longitude
ELAT top center latitude
NSTK number of point sources (subfaults) along strike
NDIP number of point sources (subfaults) down-dip
LEN segment length (km)
WID segment width (km)
STK segment strike
DIP segment dip
DTOP depth to top of fault segment (km)
SHYP along strike location (from top center) of hypocenter for this segment (km)
DHYP down-dip location (from top edge) of hypocenter for this segment (km)

**Line 4**: Number of points per fault plane

**Line 5-9**:
LON: longitude of subfault center
LAT: latitude of subfault center
DEP: depth of subfault center (km)
STK: strike
DIP: dip
AREA: area of subfault (cm2)
TINIT: initiation time when rupture reaches subfault center (sec)
DT: time step in slip velocity function (sec)
RAKE: direction of u1 axis (rake direction)
SLIP1: total slip in u1 direction (cm)
NT1: number of time points in slip rate function for u1 direction
SLIP2: total slip in u2 direction (cm)
NT2: number of time points in slip rate function for u2 direction
SLIP3: total slip in u3 (surface normal) direction (cm)
NT3: number of time points in slip rate function for u3 direction
SR1[1],…,SR1[NT1] slip rate at each time step for u1 direction (cm/sec)
SR2[1],…,SR2[NT2] slip rate at each time step for u2 direction (cm/sec)
SR3[1],…,SR3[NT3] slip rate at each time step for u3 direction (cm/sec)

Project geographic coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The geographic coordinates of source model is projected to Cartesian
coordinates wit the pre-processing tool *rconv*.

rconv -i northridge.srf -o northridge.nrf -m "+proj=merc +lon\_0=-118
+y\_0=-4050981.42 +x\_0=57329.54 +units=m +axis=enu" -x
visualization.xdmf

To find the center of fault, use *cs2cs* in *proj.4* to convert the
cooridinates:

echo -118.5150 34.3440 0.0 \| cs2cs +proj=lonlat +axis=enu +units=m +to
+proj=merc +lon\_0=-118 +axis=enu +units=m

This cooperation will project the coordinates and shift the center of
fault to the origin (0,0) in Cartesian coordinates.

Results
~~~~~~~

Source rupture starts at 7.0 s and propagates in the domain. A snapshot
of velocity is show in Figure [fig:northridge1]. The surface velocity
output is refined by subdividing each triangle into 4 subtriangles while
the domain output is not.

.. figure:: LatexFigures/1994_snap2_surface.png
   :alt: Cross-section of vertical velocity
   :width: 12.00000cm

   Cross-section of vertical velocity at surface at 7 s. The surface velocity output is refined by
   subdividing each triangle into 4 subtriangles while the domain output
   is not. The plane demonstrates the fault orientation. 




