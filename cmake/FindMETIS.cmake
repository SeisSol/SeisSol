# - Try to find METIS
# Once done this will define
#
#  METIS_FOUND          - system has METIS
#  METIS_INCLUDE_DIRS   - include directories for METIS
#  METIS_LIBRARIES      - libraries for METIS
#  METIS_64_BIT_INTEGER - METIS compiled with IDXTYPEWIDTH = 64
#
# Variables used by this module. They can change the default behaviour and
# need to be set before calling find_package:
#
#  METIS_DIR             - Where to find the base directory of metis
#  METIS_INCDIR          - Where to find the header files
#  METIS_LIBDIR          - Where to find the library files

# enforce including GKlib, if it's linked dynamically
find_package(GKlib QUIET)

find_path(METIS_INCLUDE_DIR metis.h
  HINTS ${METIS_INCLUDE_DIR} 
        ENV METIS_INCLUDE_DIR 
        ${METIS_DIR} 
        ENV METIS_DIR 
        ${PARMETIS_INCLUDE_DIR} 
        ENV PARMETIS_INCLUDE_DIR 
        ${PARMETIS_DIR} 
        ENV PARMETIS_DIR
  PATH_SUFFIXES include
  DOC "Directory where the METIS header files are located"
)

find_library(METIS_LIBRARY
  NAMES metis metis${PARMETIS_LIB_SUFFIX}
  HINTS ${PARMETIS_LIB_DIR} 
        ENV PARMETIS_LIB_DIR 
        ${PARMETIS_DIR} 
        ENV PARMETIS_DIR
        ${METIS_LIB_DIR} 
        ENV METIS_LIB_DIR 
        ${METIS_DIR} 
        ENV METIS_DIR
  PATH_SUFFIXES lib
  DOC "Directory where the METIS library is located"
)

# Try compiling and running test program
if (METIS_LIBRARY)

  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${METIS_INCLUDE_DIR})
  set(CMAKE_REQUIRED_LIBRARIES ${METIS_LIBRARY})

  # Build and run test program
  include(CheckCXXSourceRuns)
  check_cxx_source_runs("
#include <metis.h>

int main()
{
  // FIXME: Find a simple but sensible test for METIS

  return 0;
}
" METIS_TEST_RUNS)

# check whether metis was compiled with 64-bit integer
check_cxx_source_runs("
#include <metis.h>
#include <cassert>

int main()
{
  // FIXME: Find a simple but sensible test for METIS
  assert(IDXTYPEWIDTH == 64);

  return 0;
}
" _METIS_64_BIT_INTEGER)
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS REQUIRED_VARS 
                                  METIS_LIBRARY 
                                  METIS_INCLUDE_DIR 
                                  METIS_TEST_RUNS)

if(METIS_FOUND)
  set(METIS_LIBRARIES ${METIS_LIBRARY})
  set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
  if (_METIS_64_BIT_INTEGER)
    set(METIS_64_BIT_INTEGER TRUE)
  endif()
endif()

mark_as_advanced(METIS_LIBRARY METIS_INCLUDE_DIR)
