# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-io STATIC
    Datatype/Datatype.cpp
    Datatype/HDF5Type.cpp
    Datatype/MPIType.cpp
    Instance/Checkpoint/CheckpointManager.cpp
    Instance/Mesh/VtkHdf.cpp
    Instance/Metadata/Pvd.cpp
    Instance/Metadata/Xml.cpp
    Instance/Point/Csv.cpp
    Instance/Point/Pytables.cpp
    Instance/Point/TableWriter.cpp
    Reader/File/Hdf5Reader.cpp
    Reader/Distribution.cpp
    Writer/File/BinaryWriter.cpp
    Writer/File/Hdf5Writer.cpp
    Writer/Instructions/Binary.cpp
    Writer/Instructions/Data.cpp
    Writer/Instructions/Hdf5.cpp
    Writer/Instructions/Instruction.cpp
    Writer/Module/AsyncWriter.cpp
    Writer/Module/WriterModule.cpp
    Writer/Writer.cpp
    Manager.cpp
)

target_link_libraries(seissol-io PUBLIC seissol-common-properties)
