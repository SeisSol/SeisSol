..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Multiple point-sources
=======================

Standard Rupture Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~

SeisSol supports the Standard Rupture Format (SRF) for kinematic rupture
models. Details about the file format can be found at the `SCEC
Wiki <https://strike.scec.org/scecpedia/Standard_Rupture_Format>`_. With
the help of the SRF format, one may specify several subfaults (point
sources), where one may give individual source time function for each
subfault. The location, however, is given in latitude, longitude, and
depth such that it becomes necessary to convert those into the
coordinate system used by the mesh, which we call the Mesh Coordinate
System (MCS). In a software package, one usually has the option to
directly convert coordinate systems during runtime or to use an
intermediate format that does not require coordinate conversion. Here,
we opted for the second approach, because

-  we do not complicate the build and use of SeisSol and
-  we can use a binary format which greatly reduces the loading
   time.

**Hint:** Use
`krfviewer <https://github.com/SeisSol/Geodata/tree/master/krfviewer>`__
to inspect Standard Rupture Format files and to test projections.

| **Note:** There is a slight difference in SRF version 1.0 and 2.0.
| Point line 1 in 1.0: longitude latitude depth strike dip area tinit
  dt.
| Point line 1 in 2.0: longitude latitude depth strike dip area tinit dt
  vs den.

.. _netcdf-rupture-format-(nrf):

NetCDF Rupture Format (NRF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The NRF is an intermediate format for describing kinematic rupture
models. It is not meant to be used directly but it should be generated
from an SRF file. To do so, you require the tool
`rconv <https://github.com/SeisSol/SeisSol/blob/master/preprocessing/science/rconv/>`__.
Note that some python scripts required for compiling rconv are given as
a symbolic link in rconv, the link being a relative path. This means that
**you need the whole SeisSol repository to compile it**.

Specifying the MCS
^^^^^^^^^^^^^^^^^^

The main input parameter of rconv is the specification of the MCS. It is
very important to specify it right, otherwise you will get wrong moment
tensors and wrong subfault locations. The MCS can be specified by a
string describing its projection in the same way as you would use the
cartographic software `proj.4 <https://github.com/OSGeo/proj.4>`__. For
example,

::

   +proj=utm +zone=10 +datum=WGS84 +units=m +axis=ned +no_defs

would lead to UTM projection of latitude and longitude and

::

   +proj=geocent +datum=WGS84 +units=m +axis=seu +no_defs

would lead to a geocentric coordinate system. You can test the
correctness of your projection string by invoking ``cs2cs`` from the
proj.4 application suite. For example,

::

   echo 11.669 48.263 0.476 | cs2cs +proj=lonlat +datum=WGS84 +units=km +to +proj=utm +zone=33 +datum=WGS84 +units=m +axis=ned

Besides correct projections, it is also important to specify the axis
orientation with the +axis option. Examples:

-  +axis=ned: x=north, y=east, z=down
-  +axis=enu: x=east, y=north, z=up
-  +axis=seu: x=south, y=east, z=up

If the axis description does not fit your mesh, your moment tensor will
not be rotated correctly (according to strike, dip, and rake angles).


Named projections
^^^^^^^^^^^^^^^^^^

Some named projections are not recognized by proj4 (for instance, most EPSG projections). A good ressource for transposing these named projections to generic projection strings that are understood by proj4 (and rconv) can be found at
`this page <https://josm.openstreetmap.de/browser/josm/trunk/data/projection/epsg?rev=7943>`__.

Dealing with projected data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the SRF data are already projected, the projection within rconv can be by-passed using the following projection string: ``+proj=lonlat +datum=WGS84 +units=m``.


How to use Rconv
^^^^^^^^^^^^^^^^

To use e.g. Mercator projection, you can run

::

   rconv -i input.srf -o output.nrf -m “+proj=merc +lon_0=central_longtitude +y_0=translation_along_y +x_0=translation_along_x +units=m +axis=enu” -x visualization.xdmf

More projection options can be found in proj.4 website.

Checking the NRF file
^^^^^^^^^^^^^^^^^^^^^

You may manually inspect an NRF file to verify its correctness.
If you have NetCDF installed, you may enter

::

   ncdump -v subfaults sources.nrf | less

and obtain something like

::

     compound Subfault {
       double tinit ;
       double timestep ;
       double mu ;
       double area ;
       Vector3 tan1 ;
       Vector3 tan2 ;
       Vector3 normal ;
     }; // Subfault
   ...
       Subfault subfaults(source) ;
           Subfault_units subfaults:units = {"s", "s", "pascal", "m^2", "m", "m", "m"} ;
   ...
    centres = {23166.2135886125, -13374.980382819, 1250},
       {22733.1990108113, -13124.9814335668, 1250},
   ...
    subfaults =
       {10.617000000043, 0.001, 0, 250000, {-0.866025403784439, 0.5, -0}, {3.06161699786838e-17, 5.30287619362453e-17, -1}, {-0.5, -0.866025403784439, -6.12323399573677e-17}},
       {10.446600000043, 0.001, 0, 250000, {-0.866025403784439, 0.5, -0}, {3.06161699786838e-17, 5.30287619362453e-17, -1}, {-0.5, -0.866025403784439, -6.12323399573677e-17}},

which tells you the following:

-  The first point source is located at x=23166.2135886125,
   y=-13374.980382819, z=1250.
-  The STF of the source starts acting after 10.617000000043 seconds.
-  The distance between samples in the STF is 0.001 seconds.
-  The shear modulus is 0 Pa, which means that SeisSol will take the
   shear modulus from the element in which the point source resides.
-  The subfault has an area of 250000 square meters. (Be careful, in the
   SRF you have to give it in square centimetres.)
-  u_1 = {-0.866025403784439, 0.5, -0}, u_2 = {3.06161699786838e-17,
   5.30287619362453e-17, -1}, u_3 = {-0.5, -0.866025403784439,
   -6.12323399573677e-17}, where u_1 is the strike direction, u_2 is
   orthogonal to the strike direction but lies in the fault plane, and
   u_3 is the normal direction.

Using an NRF file in SeisSol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add the following section to your parameter file:

::

   &SourceType
   Type = 42
   FileName = 'sources.nrf'
   /

Pitfalls
^^^^^^^^^

Multi point-sources representation generate spurious waves at frequencies close to Vr/h with Vr the rutpure speed and h the spatial sampling of the Kinematic model.
Also, the source time function are discretized by linear interpolation, and should be adequately sampled in time to avoid sharp kinks in the source time function, which can be the source of high frequency generation.
Therefore, the kinematic model may need to be upsampled in space and/or in time, for example using this script:
https://github.com/SeisSol/SeisSol/blob/master/preprocessing/science/kinematic_models/refine_srf.py
A possible alternative is to impose the kinematic model on a dynamic rupture boundary, see :doc:`slip-rate-on-DR` for more details.
