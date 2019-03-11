Meshing with PUMGen
===================

PUMGen
------

| PUMGen
  (`https://github.com/SeisSol/PUMGen/wiki <https://github.com/SeisSol/PUMGen/wiki>`__)
  is a tool for creating and exporting meshes in an efficient format (see :ref:`PUML_mesh_format`).
  We recommand using PUMGen, in particular when large meshes have to be generated. PUMGen uses SimModeler libraries 
  for meshing.
  Here is a basic example of use:

.. code-block:: bash

   pumgen -s simmodsuite -l SimModelerLib.lic --mesh "Mesh case 1" --analysis "Analysis case 1" test.smd test

| The Analysis and meshing attributes of the model file test.smd have
  been first defined using the GUI of SimModeler, as detailed in :ref:`Meshing_with_SimModeler`.
  In particular, the --mesh and --analysis attributes are set
  according to the mesh attributes and analysis attributes defined in
  the mesh and analysis tabs (SimModeler GUI). This script will generate
  2 files: test and test.xdmf describing the mesh (see :ref:`PUML_mesh_format`).

The "Minimum insphere found" output of PUMgen allows checking the mesh
quality. The minimum insphere of a quality mesh is of
about 1/10 or more of the minimum mesh size specified. Smaller values
might impact your time step, and then your simulation time, but not so
much if using local time stepping. Also note the --analyseAR option of
PUMgen, allowing to produce an histogram of the mesh aspect ratios,
similarly to SimModeler GUI (require branch 'xml').

| PUMGen can be used to convert an ASCII \*.neu mesh file (for instance
  from gmsh) to the xmdf mesh format:

.. code-block:: bash

   pumgen test.neu test

| PUMGen can be also used to convert an netcdf mesh file (older SeisSol
  mesh file format) to the xmdf mesh format:

.. code-block:: bash

   pumgen -s netcdf test.nc test

| Multiprocessing can be used to speed up the process:

.. code-block:: bash

   mpirun -n 4 pumgen -s netcdf test.nc test

| Note that version < 10 of Simulation Modeling Suite do not support
  parallel volume meshing.

Parametrizing PUMGen with a xml file
------------------------------------

| The parametrisation of meshing and analysis attributes using the GUI
  of SimModeler can be tedious, particularly for heavy models (>Gb smd
  file, with finely sampled topography or fault surface) or when running
  parametric studies with different meshes or analysis attributes. The 'xml'
  branch of PUMGen offers a way to tag boundary conditions surfaces and 
  to set mesh attributes using a xml
  file. In addition, this allows keeping track of the meshing parameters
  in a xml file. A typical xml file can be found `here <https://github.com/TUM-I5/PUML/blob/drlts/XmlExample/meshAttributes.xml>`__.
| A typical use of the parametrisation through xml file could be:

.. code-block:: bash

   pumgen -s simmodsuite -l SimModelerLib.lic --xml meshAttributes.xml test.smd test

To determine which surface correspond to which id, the simple way of
proceeding is to use SimModeler GUI. Another option could be using the
--prbfc option of PUMgen to probe the coordinates of the model's faces.
Finally, note that with the xml option, PUMgen can also take a stl file
as an input. If used in combination with the --prbfc option of PUMGen,
it is then theorically possible to completely avoid using the GUI of
SimModeler (yet, the import from stl of SimModeler lib does not seem to
be working properly to date).
