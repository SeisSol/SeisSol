#! /usr/bin/env sh
# SPDX-FileCopyrightText: 2023-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# NOTE: just an adapted format.sh script

# sample command line usage: $0 clang-tidy $SEISSOL_SOURCE_DIR $SEISSOL_BUILD_DIR [$ARGS for clang-tidy]
# e.g. do, if SEISSOL_BUILD_DIR=$SEISSOL_SOURCE_DIR/build
# $ cd $SEISSOL_SOURCE_DIR; sh ./.ci/tidy.sh ./ ./build
# $ cd $SEISSOL_SOURCE_DIR/.ci; sh tidy.sh ../ ../build

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
        src/Model
        src/Modules
        src/Monitoring
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


    local SEISSOL_SOURCE_DIR="${1}"
    local SEISSOL_BUILD_DIR="${2}"
    shift 2

    # check for self
    if [ ! -f "${SEISSOL_SOURCE_DIR}/.ci/tidy.sh" ]; then
        echo "Please ensure that SEISSOL_SOURCE_DIR is passed as the first argument"
        exit 177
    fi

    # (we'll treat all files as headers here, but that should not have any effect on the behavior of clang-tidy)
    # regex escaping is from https://unix.stackexchange.com/a/209744
    # (it may not be 100%ly exact, but it works for our case)

    FILE_REGEX="^\$"
    for dir in ${allowlist_dir}; do
        escaped="$(printf '%s' "$dir" | sed 's/[.[\(*^$+?{|]/\\&/g')"
        FILE_REGEX="${FILE_REGEX}|${escaped}/"
    done
    for file in ${allowlist_file}; do
        escaped="$(printf '%s' "$file" | sed 's/[.[\(*^$+?{|]/\\&/g')"
        FILE_REGEX="${FILE_REGEX}|${escaped}\$"
    done

    run-clang-tidy -header-filter=$FILE_REGEX -p $SEISSOL_BUILD_DIR $@ $FILE_REGEX
}

format $@
