Preprocessing
=============

This folder contains all<sup>1</sup> tools required to set a up the input and generated source files required by SeisSol.

Tools
-----

### Meshing

- cube: This tool generates a Gambit mesh file containing a unit cube with 5n^3 tetrahedra.
  Runtime: O(n^3) = O(#elements)
  Memory usage: O(1)

### Partitioning

- cube: Generates partitioning files for the cube (see Meshing -> cube)

- gambit2seissol: This tool takes a Gambit mesh file and produces a mesh and partitioning file
  for SeisSol. It can perform one of the following or both steps:
  1. Partitioning: The mesh it partitioned into equal parts using Metis. This step is not required
     if you are running SeisSol without MPI.
  2. Reordering: The elements in the mesh are reordered. This can lead to a better performance.

<sup>1</sup>This is currently under construction. Therefore not all tools are listed here.
