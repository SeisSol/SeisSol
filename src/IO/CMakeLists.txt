add_library(seissol-io STATIC
    Datatype/Datatype.cpp
    Datatype/HDF5Type.cpp
    Datatype/MPIType.cpp
    Instance/Checkpoint/CheckpointManager.cpp
    Reader/File/Hdf5Reader.cpp
    Reader/Distribution.cpp
    Writer/File/BinaryWriter.cpp
    Writer/File/Hdf5Writer.cpp
    Writer/Instructions/Data.cpp
    Writer/Instructions/Instruction.cpp
    Writer/Module/AsyncWriter.cpp
    Writer/Module/WriterModule.cpp
    Writer/Writer.cpp
    Manager.cpp
)

target_link_libraries(seissol-io PUBLIC SeisSol-common-properties)
