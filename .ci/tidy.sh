#! /usr/bin/env sh

# NOTE: just an adapted format.sh script

# sample command line usage: $0 clang-tidy $SEISSOL_SOURCE_DIR $SEISSOL_BUILD_DIR [$ARGS for clang-tidy]
# e.g. do, if SEISSOL_BUILD_DIR=$SEISSOL_SOURCE_DIR/build
# $ cd $SEISSOL_SOURCE_DIR; sh ./.ci/tidy.sh ./ ./build
# $ cd $SEISSOL_SOURCE_DIR/.ci; sh tidy.sh ../ ../build

format() {
    # don't use a directory with whitespace
    local whitelist_dir="
        src/DynamicRupture
        "
    
    # NOTE: once the files of a directory are (almost) fully covered, consider moving it to whitelist_dir instead
    local whitelist_file="
        src/main.cpp
        "

    local SEISSOL_SOURCE_DIR="${1}"
    local SEISSOL_BUILD_DIR="${2}"
    shift 2

    # check for self
    if [ ! -f "${SEISSOL_SOURCE_DIR}/.ci/tidy.sh" ]; then
        echo "Please ensure that SEISSOL_SOURCE_DIR is passed as the first argument"
        exit 177
    fi

    # we'll treat all files as headers... Which does not matter.
    FILE_REGEX="^\$"
    for dir in ${whitelist_dir}; do
        FILE_REGEX="${FILE_REGEX}|${dir}/"
    done
    for file in ${whitelist_file}; do
        FILE_REGEX="${FILE_REGEX}|${file}\$"
    done

    python3 ${SEISSOL_SOURCE_DIR}/.ci/run-clang-tidy.py -header-filter=$FILE_REGEX -p $SEISSOL_BUILD_DIR $@ $FILE_REGEX
}

format $@
