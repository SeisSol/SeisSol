#! /usr/bin/env sh

# NOTE: This script was taken from PointCloudLibrary/pcl and adapted for SeisSol

# sample command line usage: $0 clang-format(version >= 6.0) $SEISSOL_SOURCE_DIR
# $ cd $SEISSOL_SOURCE_DIR; sh ./.ci/format.sh `which clang-format` ./
# $ cd $SEISSOL_SOURCE_DIR/.ci; sh format.sh `which clang-format` ../

format() {
    # don't use a directory with whitespace
    local allowlist_dir="
        auto-tuning/proxy/src
        src/DynamicRupture
        src/Equations/acoustic/Model
        src/Geometry
        src/Initializer/BatchRecorders
        src/Initializer/InitProcedure
        src/Initializer/Parameters
        src/Initializer/Tree
        src/IO
        src/Kernels
        src/Modules
        src/Monitoring
        src/Numerical
        src/Parallel
        src/Physics
        src/Reader
        src/ResultWriter
        src/SourceTerm
        src/tests
        "

    # NOTE: once the files of a directory are (almost) fully covered, consider moving it to allowlist_dir instead
    local allowlist_file="
        src/Equations/elastic/Model/Datastructures.h
        src/Equations/elastic/Model/IntegrationData.h
        src/Equations/viscoelastic/Model/IntegrationData.h
        src/Equations/viscoelastic2/Model/Datastructures.h
        src/Equations/viscoelastic2/Model/IntegrationData.h
        src/Equations/anisotropic/Model/Datastructures.h
        src/Equations/anisotropic/Model/IntegrationData.h
        src/Equations/poroelastic/Model/Datastructures.h
        src/Equations/poroelastic/Model/IntegrationData.h
        src/Equations/Datastructures.h
        src/Equations/Setup.h
        src/Initializer/BasicTypedefs.h
        src/Initializer/Boundary.h
        src/Initializer/DynamicRupture.h
        src/Initializer/DeviceGraph.h
        src/Initializer/GlobalData.h
        src/Initializer/GlobalData.cpp
        src/Initializer/InitialFieldProjection.h
        src/Initializer/InitialFieldProjection.cpp
        src/Initializer/InputAux.h
        src/Initializer/LTS.h
        src/Initializer/MemoryAllocator.h
        src/Initializer/MemoryAllocator.cpp
        src/Initializer/ParameterDB.h
        src/Initializer/ParameterDB.cpp
        src/Initializer/PointMapper.h
        src/Initializer/PointMapper.cpp
        src/Initializer/PreProcessorMacros.h
        src/Initializer/TimeStepping/GlobalTimestep.h
        src/Initializer/TimeStepping/GlobalTimestep.cpp
        src/Model/CommonDatastructures.h
        src/Model/Plasticity.h
        src/SeisSol.h
        src/SeisSol.cpp
        src/Main.cpp
        "


    local SEISSOL_SOURCE_DIR="${2}"
    local formatter="${1}"

    if [ ! -f "${formatter}" ]; then
        echo "Could not find clang-format. Please specify one as the first argument"
        exit 176
    fi

    local formatter_version=$(${formatter} --version)
    if [ "${formatter_version}" != "clang-format version 19.1.0" ]; then
        echo "Your clang-format tool in \"${formatter}\" does not have the correct version (should be 19.1.0). Given: ${formatter_version}"
        echo "Hint: you may install the required clang-format via pip, by typing: pip3 install clang-format==19.1.0"
        exit 176
    fi

    # check for self
    if [ ! -f "${SEISSOL_SOURCE_DIR}/.ci/format.sh" ]; then
        echo "Please ensure that SEISSOL_SOURCE_DIR is passed as the second argument"
        exit 176
    fi

    for dir in ${allowlist_dir}; do
        path=${SEISSOL_SOURCE_DIR}/${dir}
        files=$(find ${path} -type f -iname *.[ch] -o -iname *.[ch]pp -o -iname *.[ch]xx -iname *.cu)
        for file in ${files}; do
            ${formatter} -i -style=file $file
        done
    done

    for file in ${allowlist_file}; do
        ${formatter} -i -style=file ${SEISSOL_SOURCE_DIR}/$file
    done
}

format $@
