CAD models
==========

The following help pages describe how to build a structural model with
Gocad.

A simple tutorial
-----------------

See :doc:`generating-a-cad-model-using-gocad-basic-tutorial`.

Useful scripts
--------------

| `https://github.com/SeisSol/Meshing/tree/master/GocadRelatedScripts <https://github.com/SeisSol/Meshing/tree/master/GocadRelatedScripts>`__
| These python scripts are useful for creating the surfaces used in the
  CAD model. They are documented (try -h option).

Complementary information on the scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

createGOCADTSurf.py allows creating surfaces based on structured grid of
points. The lines (resp. the colums) of the point cloud should share
constant ordinates (resp. abscissa). The number of nodes on a line
(resp. on a column) should not necessary be constant, contrary to
createGOCADTSurf_NXNY.py. This script is used for creating the Sumatra
fault in the EGU Sumatra setup.

Adapting the CAD model resolution
---------------------------------

Using SimModeler
~~~~~~~~~~~~~~~~

Nowadays, high resolution topographic and bathymetric data are most of
the time available. Processing this large amount of data can really be a
challenge. For example in Gocad, intersecting such surface with other
surfaces can be time consuming and error prone. To overcome these kind
of difficulties, an idea is to coarsen the meshed topography where a
fine resolution is not necessary. For further details: :doc:`remeshing-the-topography`.
The same procedure can be also useful when the
intersection between 2 surfaces fails in gocad. In fact, remeshing one
of the surfaces can make easier the intersection in Gocad. In such a
case, all surface already intersected with the surface that we want to
remesh have to be exported to SimModeler. The mesh attributed "Use
Discrete Geometry Mesh" and "No mesh" have to be assigned to these
surfaces. This will ensure that the border nodes of the remesh surface
keep unaffected by the remeshing.

Alternative using Gocad
~~~~~~~~~~~~~~~~~~~~~~~

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
numerical precision issues in the GoCAD or in the mesher (by using
unnecessarily significant digits) For instance, for the Kaikoura
scenario, we used EPSG:3994, leading to a model centred on (6e6,-4e6) m
for a model size of roughly 500 km. It can then be a good idea to
manually centre back the model on (0,0,0).
