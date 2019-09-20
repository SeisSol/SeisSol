CAD models
==========

The following help pages describe how to build a structural model with
Gocad.

A simple tutorial
-----------------

See :doc:`generating-a-cad-model-using-gocad-basic-tutorial`.

Useful scripts
--------------

A collection of python scripts typically used to create the surfaces used in the CAD model
is available  `here <https://github.com/SeisSol/Meshing/tree/master/GocadRelatedScripts>`__.
They are documented (try -h option).
The most important script are:

-  createFaultFromCurve.py allows creating a ts surface from a fault trace. The fault trace is resampled, smoothed and extended
   using either a constant dip, a depth varying dip or an along-strike varying dip. This script has been used to generate 
   all the faults of the Kaikoura model (Ulrich et al., 2019).
-  createGOCADTSurf_NXNY.py, which allows creating a ts surface from a structured grid of points.
-  createGOCADTSurf.py, which allows creating a ts surface from a partially structured grid of points.
   Contrary to createGOCADTSurf_NXNY.py, the number of nodes on a line (resp. on a column) should not constant.
   On the other hand, the lines (resp. the colums) of the point cloud should share constant ordinates (resp. abscissa).
   This script is used for creating the Sumatra fault our Sumatra models (see Uphoff et al., 2017).
-  convertTs2Stl.py, which allows converting the geometric model from Gocad into a stl file, inputed into the mesher (e.g. SimModeler).


Processing high resolution topographic data
-------------------------------------------

High resolution topographic and bathymetric data are usually available. 
Generatign geometric models including such large dataset can be challenging.
In particular, intersecting such surface with other surfaces can be time consuming and error prone.
Here we present various strategies and tools to overcome this challenge.


Using Gdal
----------

`Gdal <https://www.gdal.org/>`__ is a powerful library to process gridded data. 
It allows for intance to easily resample or crop a dataset, and to convert file in handy file formats.
Here is a commented example of our use of gdal to create a ts surface from a high resolution topography of Nepal (file data/merged_original.tif).

.. code-block:: bash

   #resample data
   gdalwarp -s_srs EPSG:4326 -r near -tr 0.0025 0.0025 data/merged_original.tif data/file250b.tif
   #crop data
   gdalwarp -te 83.7 26. 88.1 29.4 data/file250b.tif data/file250.tif
   #change format
   gdal_translate -of netCDF -co "FOMRAT=NC4" data/file250.tif data/file250.nc
   #python script from 'GocadRelatedScript'
   #The specified hole allows to use algorithm described in 'remeshing the topography'
   python createGOCADTSurfNXNY_netcdf.py data/file250.nc data/file250.stl --proj "+init=EPSG:32645" --hole 84.8 86.5 27.1 28.3


Topographic data coarsening with SimModeler
-------------------------------------------

To avoid dealing with too large files when building the CAD model, topography data can be coarsened where
fine resolution is not necessary. For further details, see :doc:`remeshing-the-topography`.

The same procedure can be also useful when the intersection between 2 surfaces fails in gocad. In fact, remeshing one
of the surfaces can facilitate the intersection step in Gocad. In such a
case, all surface already intersected with the surface that we want to
remesh have to be exported to SimModeler. The mesh attributes "Use
Discrete Geometry Mesh" and "No mesh" have to be assigned to these
surfaces. This will ensure that the border nodes of the remesh surface
keep unaffected by the remeshing.

Alternative using Gocad
-----------------------

It can occur that the procedure described in :doc:`remeshing-the-topography`
is not applicable. For example if a first model with fine
topography has been compiled, and we want to extend it without starting
from scratch. In this case, an alternative procedure can be used:
:doc:`adapting-the-cad-model-resolution-using-gocad`.

Dealing with intersection artefacts
-----------------------------------

:doc:`manually-fixing-an-intersection-in-gocad`

On the use of projections
-------------------------

A special care must be taken when projecting from WGS84 to a projected
coordinate system (e.g. mercator) as the coordinates of the projected
model can then be centred on a point distant from (0,0), which can cause
numerical precision issues when building the geometric model or when meshing. 
For instance, for the Kaikoura
scenario, we used EPSG:3994, leading to a model centred on (6e6,-4e6) m
for a model size of roughly 500 km. It can then be a good idea to
manually centre back the model on (0,0,0).
This can usually be done by using the option +x_0=xxx and +y_0=yyy in the projection description.
