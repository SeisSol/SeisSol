Mesh Input Formats
==================

Currently, SeisSol supports two mesh formats.

.. _PUML_mesh_format:

PUML Mesh Format
~~~~~~~~~~~~~~~~

The main input format for SeisSol is the so-called PUML format,
named after the **P**arallel **U**nstructured **M**esh **L**ibrary.

The format is mostly compatible to Xdmf. For example, with a corresponding Xdmf XML file, it can be visualized in e.g. ParaView.
The data is always encoded in Hdf5. In fact, the Xdmf XML description is currently ignored by SeisSol.

The h5 file contains the data (arrays) in binary format and the xdmf describes it.
The content of the hdf5 file can be view using ``h5dump``, and the Xdmf file can be visualized using ParaView.
The boundary conditions can be extracted from a PUML mesh and visualised in ParaView using this 
`script <https://github.com/SeisSol/Meshing/blob/master/vizualizeBoundaryConditions/vizualizeBoundaryConditions.py>`_.

Here is an example of the Xdmf file (describing the hdf5 content):

.. code-block:: xml

   <?xml version="1.0" ?>
   <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
   <Xdmf Version="2.0">
    <Domain>
     <Grid Name="puml mesh" GridType="Uniform">
      <Topology TopologyType="Tetrahedron" NumberOfElements="901818">
       <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="901818 4">Sulawesi:/connect</DataItem>
      </Topology>
      <Geometry name="geo" GeometryType="XYZ" NumberOfElements="159618">
       <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="159618 3">Sulawesi:/geometry</DataItem>
      </Geometry>
      <Attribute Name="group" Center="Cell">
       <DataItem  NumberType="Int" Precision="4" Format="HDF" Dimensions="901818">Sulawesi:/group</DataItem>
      </Attribute>
      <Attribute Name="boundary" Center="Cell">
       <DataItem NumberType="Int" Precision="4" Format="HDF" Dimensions="901818">Sulawesi:/boundary</DataItem>
      </Attribute>
     </Grid>
    </Domain>
   </Xdmf>

It shows that the hdf5 file consists of 4 arrays: geometry, connect, group and boundary.

-   geometry contains the coordinates of the nodes. Dimension: (nNodes, 3)
-   connect contains the volume cell connectivity. Dimension: (nCells, 4)
-   group contains the region id of each volume cell (different properties can then be affected on different volumes). Dimension: nCells.
-   boundary contains the boundary conditions of each face of each volume cell. Dimension: nCells. (except when you use 32-bit boundary conditions)

Each tetrahedron has 4 faces. In our format, the 4 boundary condition ids (4 bits each) are store within a single integer (16 bits) variable. The values can be unpacked, for example in python using:

.. code-block:: python

   for faceId in range(0,4):
      boundaryFace[faceId] = (boundary >> (faceId*8)) & 0xFF;

| The possibles values of the boundary condition ids range from 0 to 255.
| In SeisSol, boundary conditions are historically tagged as:
| 0: regular
| 1: free surface
| 2: gravity-based free surface
| 3: dynamic rupture
| 4: dirichlet
| 5: absorbing
| 6: periodic
| 7: analytical
| n>7: dynamic rupture (See :doc:`fault-tagging`)

Here is the convention defining the face id in a tetrahedron:

.. code-block:: python

   s_vert[0,:] = [0,2,1];   s_vert[1,:] = [0,1,3];    s_vert[2,:] = [1,2,3]; s_vert[3,:] = [0,3,2];

Netcdf Input Format (Deprecated)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An older input format. It still supports periodic boundary conditions, as used e.g. for the convergence tests.
However, it requires a fixed partition encoded in the file, and is thus usually less performant than the PUML format.

Cube Generator
~~~~~~~~~~~~~~

A cube mesh generator is integrated in SeisSol as well; and it also supports periodic boundary conditions.
