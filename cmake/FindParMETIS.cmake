# SPDX-FileCopyrightText: 2023-2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-2-Clause
# SPDX-LicenseComments: Full text under /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# - Try to find ParMETIS
# Once done this will define
#
#  PARMETIS_FOUND        - system has ParMETIS
#  PARMETIS_INCLUDE_DIRS - include directories for ParMETIS
#  PARMETIS_LIBRARIES    - libraries for ParMETIS
#
# Variables used by this module. They can change the default behaviour and
# need to be set before calling find_package:
#
#  PARMETIS_DIR          - Prefix directory of the ParMETIS installation
#  PARMETIS_INCLUDE_DIR  - Include directory of the ParMETIS installation
#                          (set only if different from ${PARMETIS_DIR}/include)
#  PARMETIS_LIB_DIR      - Library directory of the ParMETIS installation
#                          (set only if different from ${PARMETIS_DIR}/lib)
#  PARMETIS_TEST_RUNS    - Skip tests building and running a test
#                          executable linked against ParMETIS libraries
#  PARMETIS_LIB_SUFFIX   - Also search for non-standard library names with the
#                          given suffix appended

#=============================================================================
# Copyright (C) 2010-2012 Garth N. Wells, Anders Logg, Johannes Ring
# and Florian Rathgeber. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# enforce including GKlib, if it's linked dynamically
find_package(GKlib QUIET)

find_path(PARMETIS_INCLUDE_DIR parmetis.h
  HINTS ${PARMETIS_INCLUDE_DIR} ENV PARMETIS_INCLUDE_DIR ${PARMETIS_DIR} ENV PARMETIS_DIR
  PATH_SUFFIXES include
  DOC "Directory where the ParMETIS header files are located"
)

find_path(METIS_INCLUDE_DIR metis.h
  HINTS ${METIS_INCLUDE_DIR} ENV METIS_INCLUDE_DIR ${METIS_DIR} ENV METIS_DIR ${PARMETIS_INCLUDE_DIR} ENV PARMETIS_INCLUDE_DIR ${PARMETIS_DIR} ENV PARMETIS_DIR
  PATH_SUFFIXES include
  DOC "Directory where the METIS header files are located"
)


find_library(PARMETIS_LIBRARY
  NAMES parmetis parmetis${PARMETIS_LIB_SUFFIX}
  HINTS ${PARMETIS_LIB_DIR} ENV PARMETIS_LIB_DIR ${PARMETIS_DIR} ENV PARMETIS_DIR
  PATH_SUFFIXES lib
  DOC "Directory where the ParMETIS library is located"
)

find_library(METIS_LIBRARY
  NAMES metis metis${PARMETIS_LIB_SUFFIX}
  HINTS ${PARMETIS_LIB_DIR} ENV PARMETIS_LIB_DIR ${PARMETIS_DIR} ENV PARMETIS_DIR
  PATH_SUFFIXES lib
  DOC "Directory where the METIS library is located"
)

# Get ParMETIS version
if(NOT PARMETIS_VERSION_STRING AND PARMETIS_INCLUDE_DIR AND EXISTS "${PARMETIS_INCLUDE_DIR}/parmetis.h")
  set(version_pattern "^#define[\t ]+PARMETIS_(MAJOR|MINOR)_VERSION[\t ]+([0-9\\.]+)$")
  file(STRINGS "${PARMETIS_INCLUDE_DIR}/parmetis.h" parmetis_version REGEX ${version_pattern})

  foreach(match ${parmetis_version})
    if(PARMETIS_VERSION_STRING)
      set(PARMETIS_VERSION_STRING "${PARMETIS_VERSION_STRING}.")
    endif()
    string(REGEX REPLACE ${version_pattern} "${PARMETIS_VERSION_STRING}\\2" PARMETIS_VERSION_STRING ${match})
    set(PARMETIS_VERSION_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
  endforeach()
  unset(parmetis_version)
  unset(version_pattern)
endif()

# Try compiling and running test program
if (PARMETIS_INCLUDE_DIR AND PARMETIS_LIBRARY AND METIS_LIBRARY)
  # Test requires MPI
  find_package(MPI QUIET REQUIRED)

  # TODO: maybe unite with FindMETIS.cmake?
  if (GKLIB_FOUND)
    set(_PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR} ${GKLIB_INCLUDE_DIRS})
    set(_PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY} ${GKLIB_LIBRARIES})
  else()
    set(_PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR})
    set(_PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
  endif()
  mark_as_advanced(_PARMETIS_INCLUDE_DIRS _PARMETIS_LIBRARIES)

  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${_PARMETIS_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
  set(CMAKE_REQUIRED_LIBRARIES ${_PARMETIS_LIBRARIES} ${MPI_LIBRARIES})

  # Build and run test program
  include(CheckCXXSourceRuns)
  check_cxx_source_runs("
#include <mpi.h>
#include <parmetis.h>

int main()
{
  // FIXME: Find a simple but sensible test for ParMETIS

  return 0;
}
" PARMETIS_TEST_RUNS)
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.2)
  find_package_handle_standard_args(ParMETIS
    REQUIRED_VARS PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR PARMETIS_TEST_RUNS
    VERSION_VAR PARMETIS_VERSION_STRING)
else()
  find_package_handle_standard_args(ParMETIS
    REQUIRED_VARS PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR PARMETIS_TEST_RUNS)
endif()

# include GKlib, if it's found
if(PARMETIS_FOUND)
  set(PARMETIS_LIBRARIES ${_PARMETIS_LIBRARIES})
  set(PARMETIS_INCLUDE_DIRS ${_PARMETIS_INCLUDE_DIRS})
endif()

mark_as_advanced(PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY METIS_LIBRARY METIS_INCLUDE_DIR)

if(PARMETIS_FOUND AND NOT TARGET ParMETIS::ParMETIS)
    add_library(ParMETIS::ParMETIS INTERFACE IMPORTED)
    set_target_properties(ParMETIS::ParMETIS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${PARMETIS_LIBRARIES}"
    )
endif()
