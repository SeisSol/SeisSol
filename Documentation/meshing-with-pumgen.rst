..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Mesh Conversion PUMGen
======================

PUMGen
------

| PUMGen
  (`https://github.com/SeisSol/PUMGen/wiki <https://github.com/SeisSol/PUMGen/wiki>`__)
  is a tool for creating and exporting meshes in an efficient format (see :ref:`PUML_mesh_format`).
  We recommend using PUMGen, in particular when large meshes have to be generated. PUMGen uses SimModeler libraries
  for meshing.
  Here is a basic example of use:

.. code-block:: bash

   pumgen -s simmodsuite -l SimModelerLib.lic --mesh "Mesh case 1" --analysis "Analysis case 1" test.smd test

| The Analysis and meshing attributes of the model file test.smd have
  been first defined using the GUI of SimModeler, as detailed in :ref:`Meshing_with_SimModeler`.
  In particular, the --mesh and --analysis attributes are set according to the mesh attributes and analysis attributes defined in the mesh and analysis tabs (SimModeler GUI).
  This script will generate 2 files: test and test.xdmf describing the mesh (see :ref:`PUML_mesh_format`).

The "Minimum insphere found" output of PUMgen allows checking the mesh
quality. The minimum insphere of a quality mesh is
about 1/10 or more of the minimum mesh size specified. Smaller values
might impact your time step, and then your simulation time, but not so
much if using local time stepping. Also note the --analyseAR option of
PUMgen, allowing to produce a histogram of the mesh aspect ratios,
similarly to SimModeler GUI.

| PUMGen can be used to convert an ASCII \*.neu mesh file (for instance
  from the GUI of SimModeler) to the xmdf mesh format:

.. code-block:: bash

   pumgen test.neu test

| PUMGen can be also used to convert a netcdf mesh file (older SeisSol
  mesh file format) to the xmdf mesh format:

.. code-block:: bash

   pumgen -s netcdf test.nc test

Or to convert a gmsh mesh file (msh2) to the xmdf mesh format:

.. code-block:: bash

   pumgen -s msh2 test.msh test

| Multiprocessing can be used to speed up the process:

.. code-block:: bash

   mpirun -n 4 pumgen -s netcdf test.nc test

| Note that version < 10 of Simulation Modeling Suite do not support
  parallel volume meshing.

Parametrizing PUMGen with an xml file
-------------------------------------

| The parametrization of meshing and analysis attributes using the GUI
  of SimModeler can be tedious, particularly for heavy models (>Gb smd
  file, with finely sampled topography or fault surface) or when running parametric studies with different meshes or analysis attributes. The --xml
  option of PUMGen offers a way to tag boundary conditions surfaces and to set mesh attributes using an xml file.
  In addition, this allows keeping track of the meshing parameters in an xml file. A typical xml file can be found `here <https://github.com/SeisSol/PUMGen/blob/master/XmlExample/meshAttributes.xml>`__.
| A typical use of the parametrization through xml file could be:

.. code-block:: bash

   pumgen -s simmodsuite -l SimModelerLib.lic --xml meshAttributes.xml test.smd test

The GUI of SimModeler can be used to find the correspondence between region/surface and id.

Velocity-aware meshing
-------------------------------------

PUMGen supports automatic mesh refinement depending on the velocity structure specified in an easi file. PUMGen generates a mesh with a local element size that satisfies the specified number of :code:`elementsPerWaveLength` for the target :code:`frequency` within the :code:`VelocityRefinementCuboid`. As a rule of thumb, running SeisSol with :code:`-DORDER=6` resolves the target frequency when using two elements per wavelength (for details see `Käser et al., 2008 <https://doi.org/10.1111/j.1365-246X.2008.03781.x>`_).

| Velocity-aware meshing is enabled within the xml file:

.. code-block:: XML

  <VelocityAwareMeshing easiFile="material.yaml" elementsPerWaveLength="2">
    <VelocityRefinementCuboid frequency="2.0" centerX="0" centerY="0" centerZ="0"
                              halfSizeX="1e5" halfSizeY="1e5" halfSizeZ="1e5"
                              bypassFindRegionAndUseGroup="1"/>
  </VelocityAwareMeshing>
