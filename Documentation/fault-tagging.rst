..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Fault tagging
=============

The dynamic rupture has the boundary tag 3, or anything larger than 64.

Using multiple fault tags enables us to initialize fault parameters segment-wise
easily. For example, when we want to model 2 segments having different dynamic friction,
we can do so by tagging them with 3 and 65. That can be then used in the fault easi file as follows:

.. code-block:: yaml

   [mu_d]: !Any
     components:
       - !GroupFilter
         groups: 3
         components: !ConstantMap
           map:
             mu_d:    0.3
       - !GroupFilter
         groups: 65
         components: !ConstantMap
           map:
             mu_d:    0.4

Currently, the only way to tag fault faces other tags than 3 with SimModeler is to use the `--xml` option of pumgen.
For example, to tag face 2 as 3 and face 8 and 9 as 65, we would
use:

.. code-block:: xml

   <boundaryCondition tag="3">2</boundaryCondition>
   <boundaryCondition tag="65">8,9</boundaryCondition>

Then pumgen is run using the xml option:

::

   pumgen -s simmodsuite -l SimModelerLib.lic --xml MeshandAnalysisAttributes.xml prefix.smd output_prefix

Note that ``<boundaryCondition tag="3">`` is equivalent to ``<dynamicRupture>``. Therefore, if you want to tag face 2 as 3, you can use either:

.. code-block:: xml

   <boundaryCondition tag="3">2</boundaryCondition>

or

.. code-block:: xml

   <dynamicRupture>2</dynamicRupture>

Note also that if a face is tagged twice, only the first tag will be considered.


Using more than 189 dynamic rupture tags
----------------------------------------

To handle more than 189 dynamic rupture tags (i.e. more than 255 boundary condition types), you will need to adjust the boundary format when building your mesh in PUMgen.

That is, add in PUMgen the option ``--boundarytype=int64`` when building your mesh.
No modification in SeisSol is needed, as it tries to infer the boundary format from the shape of the boundary array automatically.
However, to prevent mistakes with reading the format, we nevertheless recommend specifying the boundary format explicitly. (it *is* possible to confuse the boundary format, but only in some esoteric edge cases)

To do that, it suffices to specify ``pumlboundaryformat = $option`` in the ``&meshnml`` section of your SeisSol parameter file, where ``$option`` is one of the following:

- ``'auto'``: SeisSol will try to infer the boundary format automatically. This is the default option. That is, if an attribute ``boundary-format`` is present, its value will be used (``i32x4 == 0``, ``i32 == 1``, ``i64 == 2``). Otherwise, the type will be inferred from the rank of the boundary storage array: if it is two-dimensional, then we choose ``i32x4``; if it is rank 1, then we choose ``i32`` (note that ``i64`` will not be chosen automatically in this case).
- ``'i32'``: 8 bits per boundary face. That is, 189 dynamic rupture tags are possible (255 boundary condition types). It is (usually) stored as a one-dimensional 32-bit integer array (one entry per cell) in the Hdf5/binary file.
- ``'i64'``: 16 bits per boundary face. That is, 65469 dynamic rupture tags (65535 boundary condition types). It is stored as a one-dimensional 64-bit integer array (one entry per cell) in the Hdf5/binary file.
- ``'i32x4'``: 32 bits per boundary face. In short, you will have :math:`2^{32} - 65` different dynamic rupture tags available. The data is stored as a two-dimensional array (four entries per cell, one for each face) of 32-bit integers.

To try inferring which boundary format you have built your mesh for (in case it is not indicated by an attribute), you can use ``h5dump -H <yourmeshfile>.puml.h5``,
and look at the datatype and the shape of the ``boundary`` dataset. Note however, that some libraries may store ``i32`` boundaries in a 64-bit integer.
