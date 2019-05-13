.. _Meshing_with_SimModeler:

Meshing with SimModeler
=======================

The meshing workflow is presented through a simple example, by meshing
the CAD model obtained from :doc:`generating-a-cad-model-using-gocad-basic-tutorial`.
The created stl-file is imported via ``File > Import Discrete Data``.

Prerequisite
------------

The procedure to download SimModeler (GUI) and the SimModeler modeling suite (library) is detailed `here <https://github.com/SeisSol/Meshing/tree/master/SimModelerDownloadingBuilding>`__. 
Note that to be able to properly define the boundary conditions and to be able to
export the mesh in the proper format, SimModeler has to be `SeisSol
customized <https://github.com/SeisSol/Meshing/tree/master/SimModelerDownloadingBuilding#customizing-simmodeler-for-seissol>`__.

SimModeler version
------------------

We have used so far 3 main versions of SimModeler (3, 4 and 5).
Sometimes, quality meshes can be obtained on older version of SimModeler
whereas the latest version of SimModeler is not able to get quality
meshes (in that case the support of SimModeler is very reactive and
helpful). It is then important to notice that smd file created in older
versions of SimModeler can be read in all SimModeler versions. On the
other hand, smd file from the latest simModeler releases are not
backward compatible. Anyway, in most cases, we strongly recommend using
the latest version of SimModeler.

Analysis tab
------------

tab Analysis > Click twice on "New Case" on the Analysis Attributes panel.
give a name. If your SimModeler is set for SeisSol, the solver
seisSol should appear in the drop down menu.

Select the top surface (several surfaces can be selected by holding
Shift), click on the + sign > Boundary conditions > Free Surface. And
then on Apply-close (no need to enter a name).
Process similarly for the Absorbing and Dynamic rupture boundary conditions.

Meshing tab
-----------

The default Surface Meshing and Volume meshing attributes are initially
set. Their default attributes can be changed by clicking on them in the
mesh attribute tab. In particular, the Smoothing algorithm can be
changed from Laplacian to gradient (quoting SimModeler manual, "This
algorithm will generally produce better results, but at some performance
cost"). The Smoothing level can also be changed (max 1 for Volume
meshing and 4 for Surface meshing according to the manual). Finally, the
Discrete Face Rotation Angle Limit is also a parameter to consider, for
knowing to which extend the CAD model have to be matched.

| + > Mesh Size > Absolute > e.g. 5000 will define a maximum mesh size
  in the model.
| + > Gradation > Rate > e.g. 0.15 will define the coarsening rate
  within the mesh. The smaller the value, the slower the coarsening
  within the mesh.
| Click on the fault then + > Mesh Size > Absolute > e.g. 250 to define
  the on-fault size.
| Click on the fault then + > Mesh Size Propagation > propagation
  distance: e.g. 1000, scaling Factor e.g. 2. This allow the mesh to
  remain fine in a box bounding the fault. For example here, the mesh is
  coarsened away from the fault according to the gradation rate set,
  with a maximum value of 2*250m = 500m within this box.

| + > Surface Shape Metric > Aspect Ratio > e.g. 3 and
| + > Volume Shape Metric > Aspect Ratio > e.g. 6 will define quality
  level that the mesher will try to enforce. The mesher will not
  necessarily create a mesh which passes all the Shape Metric set. From
  our experience, setting additional shape metrics do not help improving
  the mesh. An easy mesh can reach AR < 10. For more complex meshes, AR
  < 40 should be expected.

Generating the mesh
-------------------

Meshing tab > Generate Mesh

Checking mesh quality
---------------------

| Display tab > Mesh Stats > check Aspect Ratio > look at extreme value
  and results spreads.
| Display tab > Region Select > change the Aspect Ratio range and
  visualize where are the badly shaped elements. If they are related to
  some geometric features of the CAD model (e.g. narrow layers, shallow
  dipping fault) then the CAD model should be modified to allow meshes
  of better quality.

Exporting the mesh
------------------

| File > Export Analysis > e.g. test.neu
| In case of a large mesh (several million cells), it is probably
  quicker to mesh using Pumgen. In this case, before meshing save the
  model file before meshing, and run the mesher using Pumgen.
