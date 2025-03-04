#! /usr/bin/env sh
# SPDX-FileCopyrightText: 2021-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# NOTE: This script was taken from PointCloudLibrary/pcl and adapted for SeisSol

# sample command line usage: $0 clang-format(version >= 6.0) $SEISSOL_SOURCE_DIR
# $ cd $SEISSOL_SOURCE_DIR; sh ./.ci/format.sh `which clang-format` ./
# $ cd $SEISSOL_SOURCE_DIR/.ci; sh format.sh `which clang-format` ../

format() {
    # don't use a directory with whitespace
    local allowlist_dir="
        src/Common
        src/DynamicRupture
        src/Equations
        src/Geometry
        src/Initializer/BatchRecorders
        src/Initializer/InitProcedure
        src/Initializer/Parameters
        src/Initializer/TimeStepping/LtsWeights
        src/IO
        src/Kernels
        src/Memory
        src/Modules
        src/Monitoring
        src/Model
        src/Numerical
        src/Parallel
        src/Physics
        src/Proxy
        src/Reader
        src/ResultWriter
        src/SourceTerm
        src/tests
        "

    # NOTE: once the files of a directory are (almost) fully covered, consider moving it to allowlist_dir instead
    local allowlist_file="
        src/Initializer/BasicTypedefs.h
        src/Initializer/CellLocalInformation.h
        src/Initializer/CellLocalMatrices.h
        src/Initializer/CellLocalMatrices.cpp
        src/Initializer/DeviceGraph.h
        src/Initializer/InitialFieldProjection.h
        src/Initializer/InitialFieldProjection.cpp
        src/Initializer/InputAux.h
        src/Initializer/ParameterDB.h
        src/Initializer/ParameterDB.cpp
        src/Initializer/PointMapper.h
        src/Initializer/PointMapper.cpp
        src/Initializer/PreProcessorMacros.h
        src/Initializer/TimeStepping/GlobalTimestep.h
        src/Initializer/TimeStepping/GlobalTimestep.cpp
        src/Solver/Estimator.h
        src/Solver/Estimator.cpp
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
