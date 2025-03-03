# SPDX-FileCopyrightText: 2023-2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

find_path(PARHIP_INCLUDE_DIR parhip_interface.h
        HINTS ${PARHIP_INCLUDE_DIR} ENV PARHIP_INCLUDE_DIR ${PARHIP_DIR} ENV PARHIP_DIR
        PATH_SUFFIXES include
        DOC "Directory where the ParHIP header files are located"
        )

find_library(PARHIP_LIBRARY
        NAMES parhip_interface parhip_interface${PARHIP_LIB_SUFFIX} parhip parhip${PARHIP_LIB_SUFFIX}
        HINTS ${PARHIP_LIB_DIR} ENV PARHIP_LIB_DIR ${PARHIP_DIR} ENV PARHIP_DIR
        PATH_SUFFIXES lib
        DOC "Directory where the ParHIP library is located"
        )

# Try compiling and running test program
if (PARHIP_INCLUDE_DIR AND PARHIP_LIBRARY)
    # Test requires MPI
    find_package(MPI QUIET REQUIRED)

    # Set flags for building test program
    set(CMAKE_REQUIRED_INCLUDES ${KAHIP_INCLUDE_DIR} ${PARHIP_INCLUDE_DIR} ${MPI_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${PARHIP_LIBRARY} ${MPI_LIBRARIES})

    # Build and run test program
    include(CheckCXXSourceRuns)
    check_cxx_source_runs("
#include <mpi.h>
#include <parhip_interface.h>
int main()
{
  // FIXME: Find a simple but sensible test for ParHIP
  return 0;
}
" PARHIP_TEST_RUNS)
endif()

# Standard package handling (for now, no versioning)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParHIP
        REQUIRED_VARS PARHIP_LIBRARY PARHIP_INCLUDE_DIR PARHIP_TEST_RUNS)

if(PARHIP_FOUND)
    set(PARHIP_LIBRARIES ${PARHIP_LIBRARY})
    set(PARHIP_INCLUDE_DIRS ${PARHIP_INCLUDE_DIR})
endif()

mark_as_advanced(PARHIP_INCLUDE_DIR PARHIP_LIBRARY)

if(PARHIP_FOUND AND NOT TARGET ParHIP::ParHIP)
    add_library(ParHIP::ParHIP INTERFACE IMPORTED)
    set_target_properties(ParHIP::ParHIP PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${PARHIP_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${PARHIP_LIBRARY}"
    )
endif()
