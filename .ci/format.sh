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
        src/Initializer/BatchRecorders
        src/Initializer/InitProcedure
        src/Initializer/Parameters
        src/Initializer/tree
        src/Modules
        src/Monitoring
        src/Numerical
        src/Parallel
        src/Physics
        src/Reader
        src/Solver
        src/SourceTerm
        src/tests
        "

    # NOTE: once the files of a directory are (almost) fully covered, consider moving it to allowlist_dir instead
    local allowlist_file="
        src/Equations/elastic/Model/datastructures.hpp
        src/Equations/elastic/Model/integrationData.hpp
        src/Equations/viscoelastic/Model/integrationData.hpp
        src/Equations/viscoelastic2/Model/datastructures.hpp
        src/Equations/viscoelastic2/Model/integrationData.hpp
        src/Equations/anisotropic/Model/datastructures.hpp
        src/Equations/anisotropic/Model/integrationData.hpp
        src/Equations/poroelastic/Model/datastructures.hpp
        src/Equations/poroelastic/Model/integrationData.hpp
        src/Equations/datastructures.hpp
        src/Equations/Setup.h
        src/Initializer/BasicTypedefs.hpp
        src/Initializer/Boundary.h
        src/Initializer/DynamicRupture.h
        src/Initializer/DeviceGraph.h
        src/Initializer/GlobalData.h
        src/Initializer/GlobalData.cpp
        src/Initializer/InitialFieldProjection.h
        src/Initializer/InitialFieldProjection.cpp
        src/Initializer/InputAux.hpp
        src/Initializer/LTS.h
        src/Initializer/MemoryAllocator.h
        src/Initializer/MemoryAllocator.cpp
        src/Initializer/ParameterDB.h
        src/Initializer/ParameterDB.cpp
        src/Initializer/PointMapper.h
        src/Initializer/PointMapper.cpp
        src/Initializer/preProcessorMacros.hpp
        src/Initializer/time_stepping/GlobalTimestep.hpp
        src/Initializer/time_stepping/GlobalTimestep.cpp
        src/Kernels/common.hpp
        src/Kernels/DynamicRupture.h
        src/Kernels/DynamicRupture.cpp
        src/Kernels/Local.h
        src/Kernels/Neighbor.h
        src/Kernels/Plasticity.h
        src/Kernels/Plasticity.cpp
        src/Kernels/PointSourceCluster.h
        src/Kernels/PointSourceClusterOnHost.h
        src/Kernels/PointSourceClusterOnHost.cpp
        src/Kernels/PointSourceClusterOnDevice.h
        src/Kernels/PointSourceClusterOnDevice.cpp
        src/Kernels/precision.hpp
        src/Kernels/Receiver.h
        src/Kernels/Receiver.cpp
        src/Kernels/Time.h
        src/Kernels/TimeCommon.h
        src/Kernels/TimeCommon.cpp
        src/Kernels/Touch.h
        src/Kernels/Touch.cpp
        src/Model/common_datastructures.hpp
        src/Model/plasticity.hpp
        src/ResultWriter/WaveFieldWriter.h
        src/ResultWriter/EnergyOutput.h
        src/ResultWriter/EnergyOutput.cpp
        src/ResultWriter/AnalysisWriter.h
        src/ResultWriter/AnalysisWriter.cpp
        src/ResultWriter/AsyncCellIDs.h
        src/ResultWriter/AsyncIO.h
        src/ResultWriter/AsyncIO.cpp
        src/ResultWriter/MiniSeisSolWriter.h
        src/ResultWriter/MiniSeisSolWriter.cpp
        src/ResultWriter/PostProcessor.h
        src/ResultWriter/PostProcessor.cpp
        src/ResultWriter/ThreadsPinningWriter.h
        src/ResultWriter/ThreadsPinningWriter.cpp
        src/SeisSol.h
        src/SeisSol.cpp
        src/main.cpp
        "


    local SEISSOL_SOURCE_DIR="${2}"
    local formatter="${1}"

    if [ ! -f "${formatter}" ]; then
        echo "Could not find clang-format. Please specify one as the first argument"
        exit 176
    fi

    local formatter_version=$(${formatter} --version)
    if [ "${formatter_version}" != "clang-format version 18.1.5" ]; then
        echo "Your clang-format tool in \"${formatter}\" does not have the correct version (should be 18.1.5). Given: ${formatter_version}"
        echo "Hint: you may install the required clang-format via pip, by typing: pip3 install clang-format==18.1.5"
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
