#! /usr/bin/env sh

# NOTE: This script was taken from PointCloudLibrary/pcl and adapted for SeisSol

# sample command line usage: $0 clang-format(version >= 6.0) $SEISSOL_SOURCE_DIR
# $ sh ./.ci/format.sh `which clang-format` ./
# $ sh format.sh `which clang-format` ../

format() {
    # don't use a directory with whitespace
    local whitelist_dir="
        src/DynamicRupture
        src/tests/DynamicRupture
        src/tests/Model
        src/tests/Reader
        src/Initializer/BatchRecorders
        src/Initializer/InitProcedure
        "
    
    # NOTE: once the files of a directory are (almost) fully covered, consider moving it to whitelist_dir instead
    local whitelist_file="
        src/Initializer/InputParameters.hpp
        src/Initializer/InputParameters.cpp
        src/Initializer/ParameterDB.h
        src/Initializer/ParameterDB.cpp
        src/Initializer/time_stepping/GlobalTimestep.hpp
        src/Initializer/time_stepping/GlobalTimestep.cpp
        src/SourceTerm/Manager.h
        src/SourceTerm/Manager.cpp
        src/SourceTerm/FSRMReader.h
        src/SourceTerm/FSRMReader.cpp
        src/Geometry/MeshReader.h
        src/Geometry/MeshReader.cpp
        src/Physics/Attenuation.hpp
        src/Physics/Attenuation.cpp
        "

    local SEISSOL_SOURCE_DIR="${2}"
    local formatter="${1}"

    if [ ! -f "${formatter}" ]; then
        echo "Could not find a clang-format. Please specify one as the first argument"
        exit 166
    fi

    # check for self
    if [ ! -f "${SEISSOL_SOURCE_DIR}/.ci/format.sh" ]; then
        echo "Please ensure that SEISSOL_SOURCE_DIR is passed as the second argument"
        exit 166
    fi

    for dir in ${whitelist_dir}; do
        path=${SEISSOL_SOURCE_DIR}/${dir}
        files=$(find ${path} -type f -iname *.[ch] -o -iname *.[ch]pp -o -iname *.[ch]xx -iname *.cu)
        for file in ${files}; do
            sed -i 's/#pragma omp/\/\/#pragma omp/g' $file
            ${formatter} -i -style=file $file
            sed -i 's/\/\/ *#pragma omp/#pragma omp/g' $file
        done
    done

    for file in ${whitelist_file}; do
        sed -i 's/#pragma omp/\/\/#pragma omp/g' $file
        ${formatter} -i -style=file $file
        sed -i 's/\/\/ *#pragma omp/#pragma omp/g' $file
    done
}

format $@
