# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

find_path(NUMA_INCLUDE_DIR numa.h
        HINTS ${NUMA_INCLUDE_DIR} ENV NUMA_INCLUDE_DIR ${NUMA_DIR} ENV NUMA_DIR
        PATH_SUFFIXES include
        DOC "Directory where the libnuma header files are located"
        )

find_library(NUMA_LIBRARY
        NAMES numa numa${NUMA_LIB_SUFFIX}
        HINTS ${NUMA_LIB_DIR} ENV NUMA_LIB_DIR ${NUMA_DIR} ENV NUMA_DIR
        PATH_SUFFIXES lib
        DOC "Directory where the libnuma is located"
        )

# Standard package handling (for now, no versioning)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NUMA
        REQUIRED_VARS NUMA_LIBRARY NUMA_INCLUDE_DIR)

if(NUMA_FOUND)
    set(NUMA_LIBRARIES ${NUMA_LIBRARY})
    set(NUMA_INCLUDE_DIRS ${NUMA_INCLUDE_DIR})
endif()

mark_as_advanced(NUMA_INCLUDE_DIR NUMA_LIBRARY)

if(NUMA_FOUND AND NOT TARGET NUMA::NUMA)
    add_library(NUMA::NUMA INTERFACE IMPORTED)
    set_target_properties(NUMA::NUMA PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${NUMA_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${NUMA_LIBRARIES}"
    )
endif()
