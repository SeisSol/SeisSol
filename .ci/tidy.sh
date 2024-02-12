#! /usr/bin/env sh

# NOTE: just an adapted format.sh script

# sample command line usage: $0 clang-tidy $SEISSOL_SOURCE_DIR $SEISSOL_BUILD_DIR [$ARGS for clang-tidy]
# e.g. do, if SEISSOL_BUILD_DIR=$SEISSOL_SOURCE_DIR/build
# $ cd $SEISSOL_SOURCE_DIR; sh ./.ci/tidy.sh ./ ./build
# $ cd $SEISSOL_SOURCE_DIR/.ci; sh tidy.sh ../ ../build

format() {
    # don't use a directory with whitespace
    local allowlist_dir="
        src/DynamicRupture
        src/Initializers/Parameters
        src/Modules
        src/Reader
        "
    
    # NOTE: once the files of a directory are (almost) fully covered, consider moving it to allowlist_dir instead
    local allowlist_file="
        src/main.cpp
        src/SeisSol.cpp
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

    python3 ${SEISSOL_SOURCE_DIR}/.ci/run-clang-tidy.py -use-color -header-filter=$FILE_REGEX -p $SEISSOL_BUILD_DIR $@ $FILE_REGEX
}

format $@
