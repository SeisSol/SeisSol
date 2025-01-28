..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

CAD models
==========

The following help pages describe how to build a structural model with either Gocad or SimModeler.

SimModeler CAD workflow
-----------------------

Since September 2019, SimModeler features powerful tools for processing discrete data (geometry in the form of meshes), which allow building structural models without the need to rely on additional CAD software.
We illustrate the SimModeler CAD workflow by building the structural model of the Palu earthquake dynamic rupture scenario (Ulrich et al., 2019).
See :doc:`simmodelerCAD-workflow`.

GOCAD CAD workflow
------------------

GOCAD is the tool we historically used to process complex geophysical data into structural models.
Since then, SimModeler developed tools for processing discrete data, in particular, a discrete surface intersection algorithm, which is much faster and more reliable than the one from GoCAD.
We, therefore, recommend the use of the SimModeler workflow. Because GoCAD may still be useful for fine processing of surface data (e.g. surface smoothing with constraints), we detail the full GoCAD workflow at: :doc:`generating-a-cad-model-using-gocad-basic-tutorial`.

Useful scripts
--------------

A collection of python scripts typically used to create the surfaces used in the CAD model
is available  `here <https://github.com/SeisSol/Meshing/tree/master/creating_geometric_models>`__.
They are documented (try -h option).
The most important scripts are:

-  ``create_fault_from_trace.py`` allows creating a meshed surface from a fault trace.
   The fault trace is resampled, smoothed, and extended using either a constant, a depth-varying, or an along-strike varying dip.
- ``create_surface_from_rectilinear_grid.py`` allows creating a meshed surface from a (possibly sparse, e.g. Slab2.0 dataset) structured dataset (e.g. netcdf).
-  ``create_surface_from_structured_grid.py`` allows creating a meshed surface from structured grid of nodes.
   Vertices in a row or a column do not necessary share the same x and y values (non rectilinear grid).
-  ``surface_from_one_dim_structured_grid.py`` allows creating a meshed surface from a partially structured grid of nodes.
   The point set should consist of serveral lines of nodes of same y (or x, depending on args.axis) coordinates.
   Contrary to ``create_surface_from_structured_grid.py``, the number of nodes on a line (resp. on a column) is not necessarily constant.
   On the other hand, the lines (resp. the columns) of the point cloud should share constant ordinates (resp. abscissa).
-  ``convertTs.py`` allows converting the geometric model from Gocad into another supported format (e.g. stl, bstl).


Processing high-resolution topographic data
-------------------------------------------

High resolution topographic and bathymetric data are usually available.
Generating geometric models including such large datasets can be challenging.
In particular, intersecting such surfaces with other surfaces can be time-consuming and error-prone.
Here we present various strategies and tools to overcome this challenge.


Using Gdal
----------

`Gdal <https://www.gdal.org/>`__ is a powerful library to process gridded data.
It allows, for instance, to easily resample or crop a dataset, and to convert files in handy file formats.
Here is a commented example of our use of Gdal to create a ts surface from a high-resolution topography of Nepal (file data/merged_original.tif).

.. code-block:: bash

   #resample data
   gdalwarp -s_srs EPSG:4326 -r near -tr 0.0025 0.0025 data/merged_original.tif data/file250b.tif
   #crop data
   gdalwarp -te 83.7 26. 88.1 29.4 data/file250b.tif data/file250.tif
   #change format
   gdal_translate -of netCDF -co "FOMRAT=NC4" data/file250.tif data/file250.nc
   #python script from 'creating_geometric_model'
   #The specified hole allows to use algorithm described in 'remeshing the topography'
   python3 create_surface_from_rectilinear_grid.py data/file250.nc data/file250.stl --proj "+init=EPSG:32645" --hole 84.8 86.5 27.1 28.3


Topographic data coarsening with SimModeler
-------------------------------------------

To avoid dealing with too large files when building the CAD model, topography data can be coarsened where
fine resolution is not necessary. For further details, see :doc:`remeshing-the-topography`.

The same procedure can be also useful when the intersection between 2 surfaces fails in Gocad.
In fact, creating a clean mesh of one of the surfaces can facilitate the intersection step in Gocad. In such a
case, all surface already intersected with the surface that we want to
mesh again have to be exported to SimModeler. The mesh attributes "Use
Discrete Geometry Mesh" and "No mesh" have to be assigned to these
surfaces. This will ensure that the border nodes of the new meshed surfaces keep unchanged.

Alternative using Gocad
-----------------------

It can occur that the procedure described in :doc:`remeshing-the-topography`
is not applicable. For example, if a first model with fine
topography has been compiled, and we want to extend it without starting
from scratch. In this case, an alternative procedure can be used:
:doc:`adapting-the-cad-model-resolution-using-gocad`.

Dealing with intersection artifacts
-----------------------------------

:doc:`manually-fixing-an-intersection-in-gocad`

.. _On the use of projections:

On the use of projections
-------------------------

Special care must be taken when projecting from WGS84 to a projected
coordinate system (e.g. Mercator) as the coordinates of the projected
model can then be centered on a point distant from (0,0), which can cause
numerical precision issues when building the geometric model or when meshing.
For instance, for the Kaikoura
scenario, we used EPSG:3994, leading to a model centered on (6e6,-4e6) m
for a model size of roughly 500 km. It can then be a good idea to
manually center back the model on (0,0,0).
This can usually be done by using the option +x_0=xxx and +y_0=yyy in the projection description.
