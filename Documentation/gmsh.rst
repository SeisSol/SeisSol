Gmsh
====

Introduction
------------

`Gmsh <http://gmsh.info/>`_ is an open-source
mesher, able to generate quality meshes for not too complex models.
Contrary to SimModeler, it can build geometric models from scratch. It
is particularly useful for simple idealised models, e.g. planar fault,
no topography. Two examples of model generation using Gmsh are provided
at this `link <https://github.com/SeisSol/SeisSol/tree/master/preprocessing/meshing/gmsh_example>`_.
The purpose of this tutorial is not to explain all functions of gmsh,
but to describe the features useful for setting a model and getting a
mesh for SeisSol.

Coarsening the mesh
-------------------

The Attractor Field feature is used to coarsen the mesh away from the
fault. As it requires a "Ruled Surface", we define an artificial "Ruled
Surface" used only for that feature. The rate of coarsening is defined
empirically by combining a linear increase of the mesh size in the near
field and a quadratic increase in the far size. Example:

.. code-block:: YAML

  Field[2].F = Sprintf("0.1*F1 +(F1/5.0e3)^2 + %g", lc_fault);

Boundary conditions
-------------------

The free-surface (resp. dynamic rupture, resp. absorbing) boundary
conditions are set using Physical Surfaces 101 (resp. 103, resp. 105). 
The volumes should also be put into Physical Volume to be exported into
the mesh.Here is an example from tpv33:

.. code-block:: YAML

  Physical Surface(101) = {252, 258, 260, 262};
  Physical Surface(105) = {242, 244, 246, 248, 250, 254, 256};
  Physical Surface(103) = {2, 3, 4};
  Physical Volume(2) = {276,278};
  Physical Volume(3) = {274};

Generating the mesh
-------------------

| Once the geometry and boundary conditions are set, the mesh can be
  obtained using the following command:

.. code-block:: YAML

  gmsh test.geo -3 -optimize -format neu

Note that the '-format neu' is only possible since gmsh 4.0
For previous versions, we used `gmsh2gambit <https://github.com/SeisSol/SeisSol/tree/master/preprocessing/meshing/gmsh2gambit>`_
on the msh mesh generated with:

.. code-block:: YAML

  gmsh test.geo -3 -optimize

| The optimize option is very important. If not used, mesh of very poor
  quality are generated. Optimizing the mesh using "gmsh optimize", then
  "net_gen optimize" and finally "gmsh optimize" may allow getting
  meshes of higher quality.

gmsh to SimModeler
------------------

| It is possible to create the geometry with gmsh and then mesh it with
  SimModeler. A way of doing so it to put all surfaces of the model in a
  "physical surface", mesh them (-2) and output them to a stl file (e.g.
  -o test.stl). Then the stl file can be opened with SimModeler and the
  mesh can be generated.
| If SimModeler merges some features of the geometry, it is then
  necessary to isolate the features in different stl files (i.e. running
  several times ``gmsh ___.geo -2 -o ___.stl`` with different surfaces
  listed in the physical surface listing). Then the solid name attribute
  of the stl files have to be modified. Finally, the stl files can be
  merged into a single stl file, to be opened in SimModeler.

mirroring a mesh
----------------

In order to get a maximum accuracy, it is sometimes necessary (e.g. for
benchmarks) to mirror a mesh. To get a mirrored mesh, a half mesh is
first generated. Then it is converted to netcdf format (one partition),
using PUMgen. Finally, this
`matlab script <https://github.com/SeisSol/SeisSol/blob/master/preprocessing/meshing/mirror_mesh.m>`_
allows creating the mirrored mesh.
