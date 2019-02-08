meshing partionning with PUMgen (deprecated)
============================================

The material on this page is deprecated: SeisSol now uses an new
unpartionned mesh format, partitionning being done during simulation.
Yet, optimized partitionning for local time stepping is not implemented
on the 'master' branch, but only on the 'easi' branch, which is our
current developpement branch. As a consequence, the meshing/partionning
workflow with PUMgen can still be of interest and is therefore on this
wiki.

Introduction
------------

| PUMgen
  (`https://github.com/TUM-I5/PUML/wiki <https://github.com/TUM-I5/PUML/wiki>`__)
  is a tool written by S. Rettenberger allowing to mesh, partition and
  export the mesh file in an efficient netcdf format. We recommand using
  PUMgen, in particular for large meshes. The meshing is done through
  SimModeler libraries and the partition through Metis.
| Here is a meshing and partitioning example using PUMgen:
| ``pumgen -s simmodsuite -l SimModelerLib.lic --mesh "Mesh case 1" --analysis "Analysis case 1" test.smd 20 test.20.nc``
| The Analysis and meshing attributes of the model file test.smd have
  been first defined using the GUI of SimModeler, as previously
  presented. In particular, the --mesh and --analysis attributes are set
  according to the mesh attributes and analysis attributes defined in
  the mesh and analysis tabs (SimModeler GUI). 20 is the number of
  partition required. Note it can be convenient to also write the number
  of partitions in the file name.
| The "Minimum insphere found" output by PUMgen allow you checking the
  quality of the mesh. The minimum insphere of a quality mesh should be
  of about 1/10 or more of the minimum mesh size specified. Smaller
  values might impact your time step, and then your simulation time, but
  not so much if using local time stepping. Also note the
  ``--analyseAR`` option of PUMgen, allowing to produce an histogram of
  the mesh aspect ratios, similarly to SimModeler GUI.
| PUMgen can be used to partition a \*.neu file:
| ``pumgen test.neu 2 test.2.nc``
| or to repartion a netcdf file:
| ``pumgen -s netcdf test.4.nc 20 test.20.nc``
| Mpirun can be use to speed up the process:
| ``mpirun -n 4 pumgen -s netcdf planar.4.nc 20 planar.20.nc``
| Note that version < 10 of Simulation Modeling Suite do not support
  parallel volume meshing.

Parametrizing PUMgen with a xml file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| The parametrisation of meshing and analysis attributes using the GUI
  of SimModeler can be tedious, particularly for heavy models (>Gb smd
  file, with finely sampled topography or fault surface) or when running
  parametric studies with different meshes or mesh attributes. PUMgen
  offers a way to tag surfaces with the correct boundary condition, and
  to set the various mesh attributes using a xml file. In addition, this
  allows keeping track of the meshing parameters in a small ASCII file.
  A typical xml file can be found here:
| `https://github.com/TUM-I5/PUML/blob/drlts/XmlExample/meshAttributes.xml <https://github.com/TUM-I5/PUML/blob/drlts/XmlExample/meshAttributes.xml>`__
| A typical use of the parametrisation through xml file could be:
| ``pumgen -s simmodsuite -l SimModelerLib.lic --xml  meshAttributes.xml test.smd 2 test.2.nc``
| To determine which surface correspond to which id, the simple way of
  proceeding is to use SimModeler GUI. Another option could be using the
  ``--prbfc`` option of PUMgen to probe the coordinates of the model's
  faces. Here is an example of output:

   | Thu Jul 13 16:33:15, Info: There are 2 regions in the model
   | Thu Jul 13 16:33:15, Info: There are 8 faces on model region 1 :
   | Thu Jul 13 16:33:15, Info: 1
   | Thu Jul 13 16:33:15, Info: 4
   | Thu Jul 13 16:33:15, Info: 8
   | Thu Jul 13 16:33:15, Info: 2
   | Thu Jul 13 16:33:15, Info: 11
   | Thu Jul 13 16:33:15, Info: 3
   | Thu Jul 13 16:33:15, Info: 5
   | Thu Jul 13 16:33:15, Info: 6
   | Thu Jul 13 16:33:15, Info: There are 5 faces on model region 2 :
   | Thu Jul 13 16:33:15, Info: 6
   | Thu Jul 13 16:33:15, Info: 11
   | Thu Jul 13 16:33:15, Info: 9
   | Thu Jul 13 16:33:15, Info: 10
   | Thu Jul 13 16:33:15, Info: 7
   | Thu Jul 13 16:33:15, Info: Face information:
   | Thu Jul 13 16:33:15, Info: There are 1403 polygons and 959 points
     on model face 1 , e.g.:
   | Thu Jul 13 16:33:15, Info: Polygon 0 has the following points: 0 1
     2
   | Thu Jul 13 16:33:15, Info: Point 0 : ( 6.14593e+06 , -3.9368e+06 ,
     -20000 )
   | Thu Jul 13 16:33:15, Info: Point 1 : ( 6.14624e+06 , -3.93588e+06 ,
     -20000 )
   | Thu Jul 13 16:33:15, Info: Point 2 : ( 6.14719e+06 , -3.93609e+06 ,
     -18893.5 )
   | Thu Jul 13 16:33:15, Info: There are 10148 polygons and 6241 points
     on model face 2 , e.g.:
   | Thu Jul 13 16:33:15, Info: Polygon 0 has the following points: 0 1
     2
   | Thu Jul 13 16:33:15, Info: Point 0 : ( 6.20001e+06 , -3.86484e+06 ,
     -22000.1 )
   | Thu Jul 13 16:33:15, Info: Point 1 : ( 6.19832e+06 , -3.86431e+06 ,
     -22970.5 )
   | Thu Jul 13 16:33:15, Info: Point 2 : ( 6.19901e+06 , -3.86362e+06 ,
     -22966.8 )

Finally, note that with the xml option, PUMgen can also take a stl file
as an input. If used in combination with the ``--prbfc`` option of
PUMgen, it is then theorically possible to completely avoid using the
GUI of SimModeler (yet, the import from stl of SimModeler lib does not
seem to be working properly to date).

Partionning with PUMgen for LTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| When using Local time stepping, the partition has to take into account
  that there are different time steps, i.e. small elements are more
  expensive than big elements. When using ClusteredLTS=2, an adapted
  partition is obtained using the -v 2 option of PUMgen (compiled from
  the lts_weights branch). Example:
| ``pumgen -v 2 test.neu 2 test.2.nc``
| The dynamic rupture cells are more expensive than the normal cells.
  This can be accounted for using option '--dr-to-cell-ratio 1'.
| Finally, the cell time step depends on its material properties. Then,
  if the medium is heterogeneous, the partitioning can me improved by
  taking the velocity model into account. The `following
  commit <https://github.com/TUM-I5/PUML/commit/ecf51964eb81ee7d721bbb2e89b88f9a85493104>`__
  gives a good example for adding a custom velocity model into PUMgen.
  The velocity-model can then be invoked in PUMgen using the
  velocity-model option:
| ``pumgen --dr-to-cell-ratio 1 -v 2 --velocity-model sumatra1223 -s netcdf test.32.nc 32 newmesh.32.nc``
