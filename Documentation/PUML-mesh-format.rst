..
  SPDX-FileCopyrightText: 2019-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Mesh Input Formats
==================

Currently, SeisSol supports two mesh formats.

.. _PUML_mesh_format:

PUML Mesh Format
~~~~~~~~~~~~~~~~

The main input format for SeisSol is called PUML;
an acronym for **P**\ arallel **U**\ nstructured **M**\ esh **L**\ ibrary.

The data is stored in a single Hdf5 file, consisting of the following datasets:
-   ``/geometry``: contains the geometric coordinates of the vertices. Dimension: (nNodes, 3)
-   ``/connect``: describes the tetrahedra by their vertices. Dimension: (nCells, 4)
-   ``/group``: contains the group ID for each tetrahedron; which is used in selecting different material properties. Dimension: (nCells,).
-   ``/boundary``: contains the boundary condition information for all four faces of each tetrahedron. Dimension: (nCells,) or (nCells, 4).

The content of the Hdf5 file can be view using ``h5dump``; most notably ``h5dump -H <file>`` will output only the dataset headers.

The format itself is mostly compatible to Xdmf, given a corresponding Xdmf XML file. With such, it can be visualized in e.g. ParaView.
Other than that, the Xdmf XML description is currently ignored by SeisSol.

The boundary conditions can be extracted from a PUML mesh and visualised in ParaView using this
`script <https://github.com/SeisSol/Meshing/blob/master/vizualizeBoundaryConditions/vizualizeBoundaryConditions.py>`_.

An example Hdf5 file looks as follows:

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

It shows that the hdf5 file consists of the 4 arrays: geometry, connect, group and boundary.

In the default format (``int32``), the 4 boundary condition ids for each tetrahedron (8 bits each) are store within a single integer (32 bits) variable. The values can be unpacked, for example in python using:

.. code-block:: python

   for faceId in range(0,4):
      boundaryFace[faceId] = (boundary >> (faceId*8)) & 0xFF;

Other boundary formats (``int64``) have a 16-bit offset and use 0xffff as a mask instead. The format ``int32x4`` stores each boundary value in an array value of its own, instead of compressing all four values into one integer.

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
| n>64: dynamic rupture (See :doc:`fault-tagging`)

The following convention for defining a face ID is used:

.. code-block:: python

   s_vert[0,:] = [0,2,1];   s_vert[1,:] = [0,1,3];    s_vert[2,:] = [1,2,3]; s_vert[3,:] = [0,3,2];

Netcdf Input Format (Deprecated)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An older input format. It still supports periodic boundary conditions, as used e.g. for the convergence tests.
However, it requires a fixed partition encoded in the file, and is thus usually far more restricted and less performant than the PUML format.

Cube Generator
~~~~~~~~~~~~~~

A cube mesh generator is integrated in SeisSol as well; it also supports periodic boundary conditions.
