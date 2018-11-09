Meshing with PUMGen
===================

Mesh format
-----------

SeisSol has recently gained in flexibility and user-friendliness through
a new mesh format. Previously, SeisSol used a netcdf mesh format which
contained the partitioning of the mesh, which was computed using a
preprocessing algorithm. This netcf format was then partition-dependent.
In addition, the velocity model was required when partitioning for local
time stepping, which was source of complexity. As a consequence, we have
now switched to a non-partioned xdmf mesh format, which is more flexible
and lighter. Partitioning is now done during simulation. Note that an
optimized partitioning for local time stepping is not yet available on
the 'master' branch, but only on the 'easi' branch, which is our current
development branch (that we strongly encourage to use). The workflow
allowing generating partitionned mesh in the older netcdf format is
described here: [[meshing partionning with PUMgen (deprecated)]]. Note
that netcdf can be easily converted to the new mesh format using PUMGen.

PUMGen
------

| PUMGen
  (`https://github.com/SeisSol/PUMGen/wiki <https://github.com/SeisSol/PUMGen/wiki>`__)
  is a tool for creating and exporting meshes in an efficient xdmf file.
  We recommand using PUMGen, in particular for large meshes. The meshing
  is done through SimModeler libraries and the partition through Metis.
  Here is a basic example of use:
| ``pumgen -s simmodsuite -l SimModelerLib.lic --mesh "Mesh case 1" --analysis "Analysis case 1" test.smd test``
| The Analysis and meshing attributes of the model file test.smd have
  been first defined using the GUI of SimModeler, as previously
  presented. In particular, the --mesh and --analysis attributes are set
  according to the mesh attributes and analysis attributes defined in
  the mesh and analysis tabs (SimModeler GUI). This script will generate
  2 files: test and test.xdmf describing the mesh.

The "Minimum insphere found" output by PUMgen allow you checking the
quality of the mesh. The minimum insphere of a quality mesh should be of
about 1/10 or more of the minimum mesh size specified. Smaller values
might impact your time step, and then your simulation time, but not so
much if using local time stepping. Also note the --analyseAR option of
PUMgen, allowing to produce an histogram of the mesh aspect ratios,
similarly to SimModeler GUI (require branch 'xml').

| PUMGen can be used to convert an ASCII \*.neu mesh file (for instance
  from gmsh) to the xmdf mesh format:
| ``pumgen test.neu test``
| PUMGen can be also used to convert an netcdf mesh file (older SeisSol
  mesh file format) to the xmdf mesh format:
| ``pumgen -s netcdf test.nc test``
| Mpirun can be use to speed up the process:
| ``mpirun -n 4 pumgen -s netcdf test.nc test``
| Note that version < 10 of Simulation Modeling Suite do not support
  parallel volume meshing.

Parametrizing PUMGen with a xml file
------------------------------------

| The parametrisation of meshing and analysis attributes using the GUI
  of SimModeler can be tedious, particularly for heavy models (>Gb smd
  file, with finely sampled topography or fault surface) or when running
  parametric studies with different meshes or mesh attributes. The 'xml'
  branch of PUMGen offers a way to tag surfaces with the correct
  boundary condition, and to set the various mesh attributes using a xml
  file. In addition, this allows keeping track of the meshing parameters
  in a small ASCII file. A typical xml file can be found here:
| `https://github.com/TUM-I5/PUML/blob/drlts/XmlExample/meshAttributes.xml <https://github.com/TUM-I5/PUML/blob/drlts/XmlExample/meshAttributes.xml>`__
| A typical use of the parametrisation through xml file could be:
| ``pumgen -s simmodsuite -l SimModelerLib.lic --xml meshAttributes.xml test.smd test``

To determine which surface correspond to which id, the simple way of
proceeding is to use SimModeler GUI. Another option could be using the
--prbfc option of PUMgen to probe the coordinates of the model's faces.
Finally, note that with the xml option, PUMgen can also take a stl file
as an input. If used in combination with the --prbfc option of PUMGen,
it is then theorically possible to completely avoid using the GUI of
SimModeler (yet, the import from stl of SimModeler lib does not seem to
be working properly to date).
