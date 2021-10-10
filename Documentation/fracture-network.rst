.. _SeisSol-fracture-network:

SeisSol for Fracture/Fault Network
==================================

Introduction
------------

SeisSol can handle complex geometry (Wolherr et al. 2018; Wolherr et al. 2019; Ulrich et al. 2019). It is possible to use branched faults, dipping fault geometry, and crosscutting fault (Palgunadi et al. 2020). It is possible to use branched faults, dipping fault geometry, and crosscutting fault (Palgunadi et al. 2020). However, to work on fracture/fault network simulation with individual face tagging (see `fault tagging <https://seissol.readthedocs.io/en/latest/fault-tagging.html>`__), the `SeisSol-master-branch <https://github.com/SeisSol/SeisSol>`__ cannot handle for number of fracture/fault segments more than 189 (see issue `270 <https://github.com/SeisSol/SeisSol/issues/274>`__). Therefore, several modifications in pumgen (https://seissol.readthedocs.io/en/latest/meshing-with-pumgen.html), PUML submodule (https://seissol.readthedocs.io/en/latest/PUML-mesh-format.html), and SeisSol should be performed.

Modification in PUMGen
----------------------

The current `PUMGen-master-branch <https://github.com/SeisSol/PUMGen>`__ can handle only for fracture/fault tagging less than 189. For more fracture numbers, one should clone from `PUMGen64bit <https://github.com/palgunadi1993/PUMGen/tree/PUMGen64bit>`__ (https://github.com/palgunadi1993/PUMGen.git) or compare the necessary changes in `changes <https://github.com/palgunadi1993/PUMGen/commit/75e7967f593049f4638e9790be4b676c092a9662>` The compilation still follows the installation tutorial (https://github.com/SeisSol/PUMGen/wiki/How-to-compile-PUMGen).

As the only way to tag fault/fracture faces is to use the -xml option of pumgen, we should write:

.. code-block:: xml

   <boundaryCondition tag="3">13245</boundaryCondition>
   .
   .
   .
   <boundaryCondition tag="900">12345,14325</boundaryCondition>

The pumgen is run using xml option (https://seissol.readthedocs.io/en/latest/fault-tagging.html).

Modification in PUML
--------------------

As we change the mesh data type from 32bit to 64bit, we must change the mesh reader (for PUML mesh format using Simmodeler software) in SeisSol. The PUML mesh reader is located in submodules/PUML/PUML.h. Therefore, the necessary changes can be viewed from `PUML64bit <https://github.com/palgunadi1993/PUML2/commit/392115f10d1aa774865dd927e50a8a9bfbdf5ed1>`__. 

Modification in SeisSol
-----------------------

There are several changes needed in SeisSol in accordance to the changes in mesh data type. The list of files are as the following:
-  src/Geometry/PUMLReader.cpp
-  src/Initializer/ParameterDB.cpp
-  src/Initializer/time_stepping/LtsWeights.cpp
-  src/Initializer/time_setpping/LtsWeights.h

The necessary changes can be found in `SeisSol64bit <https://github.com/palgunadi1993/SeisSol/commit/e9e4cb65ca86460fdbb5636cd0fabfaec5221968>`__. One might want to clone from `SeisSol64bit-branch <https://github.com/palgunadi1993/SeisSol/tree/SeisSol64bit>`__ (https://github.com/palgunadi1993/SeisSol.git).
