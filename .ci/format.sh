#! /usr/bin/env sh

# NOTE: This script was taken from PointCloudLibrary/pcl and adapted for SeisSol

# sample command line usage: $0 clang-format(version >= 6.0) $SEISSOL_SOURCE_DIR
# $ cd $SEISSOL_SOURCE_DIR; sh ./.ci/format.sh `which clang-format` ./
# $ cd $SEISSOL_SOURCE_DIR/.ci; sh format.sh `which clang-format` ../

format() {
    # don't use a directory with whitespace
    local allowlist_dir="
        src/DynamicRupture
        src/Geometry
        src/Initializer/Parameters
        src/tests/Common
        src/tests/DynamicRupture
        src/tests/Initializer
        src/tests/Kernel
        src/tests/Model
        src/tests/Reader
        src/tests/SourceTerm
        src/Initializer/BatchRecorders
        src/Initializer/InitProcedure
        src/Modules
        src/Monitoring
        src/Physics
        src/SourceTerm
        src/Reader
        "
    
    # NOTE: once the files of a directory are (almost) fully covered, consider moving it to allowlist_dir instead
    local allowlist_file="
        src/Equations/elastic/Model/datastructures.hpp
        src/Equations/viscoelastic2/Model/datastructures.hpp
        src/Equations/anisotropic/Model/datastructures.hpp
        src/Equations/poroelastic/Model/datastructures.hpp
        src/Initializer/BasicTypedefs.hpp
        src/Initializer/InputAux.hpp
        src/Initializer/ParameterDB.h
        src/Initializer/ParameterDB.cpp
        src/Initializer/preProcessorMacros.hpp
        src/Initializer/time_stepping/GlobalTimestep.hpp
        src/Initializer/time_stepping/GlobalTimestep.cpp
        src/Initializer/tree/LTSSync.hpp
        src/Kernels/common.hpp
        src/Kernels/PointSourceCluster.h
        src/Kernels/PointSourceClusterOnHost.h
        src/Kernels/PointSourceClusterOnHost.cpp
        src/Kernels/PointSourceClusterOnDevice.h
        src/Kernels/PointSourceClusterOnDevice.cpp
        src/Kernels/Touch.h
        src/Kernels/Touch.cpp
        src/Model/common_datastructures.hpp
        src/Model/plasticity.hpp
        src/Parallel/AcceleratorDevice.h
        src/Parallel/AcceleratorDevice.cpp
        src/Parallel/DataCollector.h
        src/Parallel/Helper.hpp
        src/ResultWriter/WaveFieldWriter.h
        src/ResultWriter/EnergyOutput.h
        src/ResultWriter/EnergyOutput.cpp
        src/SeisSol.h
        src/SeisSol.cpp
        src/main.cpp
        "

    local SEISSOL_SOURCE_DIR="${2}"
    local formatter="${1}"

    if [ ! -f "${formatter}" ]; then
        echo "Could not find a clang-format. Please specify one as the first argument"
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
