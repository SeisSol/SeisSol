# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-2-Clause

# Adapted from FindMETIS.cmake which in turn seems to be adapted from FindParMETIS.cmake
# Hence, BSD-2-Clause (not BSD-3-Clause)

# - Try to find GKLIB
# Once done this will define
#
#  GKLIB_FOUND          - system has GKLIB
#  GKLIB_INCLUDE_DIRS   - include directories for GKLIB
#  GKLIB_LIBRARIES      - libraries for GKLIB
#  GKLIB_64_BIT_INTEGER - GKLIB compiled with IDXTYPEWIDTH = 64
#
# Variables used by this module. They can change the default behaviour and
# need to be set before calling find_package:
#
#  GKLIB_DIR             - Where to find the base directory of gklib
#  GKLIB_INCDIR          - Where to find the header files
#  GKLIB_LIBDIR          - Where to find the library files

# we take one of the header files here and check for it

find_path(GKLIB_INCLUDE_DIR gk_defs.h
  HINTS ${GKLIB_INCLUDE_DIR} 
        ENV GKLIB_INCLUDE_DIR
        ${GKLIB_DIR} 
        ENV GKLIB_DIR
  PATH_SUFFIXES include
  DOC "Directory where the GKlib header files are located"
)

find_library(GKLIB_LIBRARY
  NAMES GKlib
  HINTS ${GKLIB_LIB_DIR} 
        ENV GKLIB_LIB_DIR
        ${GKLIB_DIR} 
        ENV GKLIB_DIR
  PATH_SUFFIXES lib
  DOC "Directory where the GKlib library is located"
)

# Try compiling and running test program
if (GKLIB_LIBRARY)
  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${GKLIB_INCLUDE_DIR})
  set(CMAKE_REQUIRED_LIBRARIES ${GKLIB_LIBRARY})

  # Build and run test program
  include(CheckCXXSourceRuns)
  check_cxx_source_runs("
#include <gk_defs.h>

int main()
{
  // FIXME: Find a simple but sensible test for GKlib

  return 0;
}
" GKLIB_TEST_RUNS)
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GKlib REQUIRED_VARS 
                                  GKLIB_LIBRARY 
                                  GKLIB_INCLUDE_DIR 
                                  GKLIB_TEST_RUNS)

if(GKLIB_FOUND)
  set(GKLIB_LIBRARIES ${GKLIB_LIBRARY})
  set(GKLIB_INCLUDE_DIRS ${GKLIB_INCLUDE_DIR})
endif()

mark_as_advanced(GKLIB_LIBRARY GKLIB_INCLUDE_DIR)
