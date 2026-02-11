# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

function(capitalize INARG OUTARG)
  string(REPLACE "_" ";" PARTS ${INARG})
  set(PREOUTARG "")
  foreach(PART IN LISTS PARTS)
    string(SUBSTRING ${PART} 0 1 FIRST)
    string(TOUPPER ${FIRST} FIRSTU)
    string(LENGTH ${PART} LEN)
    string(SUBSTRING ${PART} 1 ${LEN} NEXT)
    string(PREPEND NEXT ${FIRSTU})
    string(APPEND PREOUTARG ${NEXT})
  endforeach()
  set(${OUTARG} ${PREOUTARG} PARENT_SCOPE)
endfunction()

set (BUILDTYPE_STR "Cpu")
set (DEVICE_BACKEND_STR "None")
if (WITH_GPU)
  set (BUILDTYPE_STR "Gpu")
  if (DEVICE_BACKEND STREQUAL "cuda")
    set (DEVICE_BACKEND_STR "Cuda")
  endif()
  if (DEVICE_BACKEND STREQUAL "hip")
    set (DEVICE_BACKEND_STR "Hip")
  endif()
  if (DEVICE_BACKEND STREQUAL "hipsycl" OR DEVICE_BACKEND STREQUAL "acpp" OR DEVICE_BACKEND STREQUAL "oneapi")
    set (DEVICE_BACKEND_STR "Sycl")
  endif()
endif()

configure_file("Alignment.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/Alignment.h")

if (EQUATIONS STREQUAL "viscoelastic")
  set(PARAMETER_MATERIAL "Viscoelastic")
  set(PARAMETER_VISCOMODE "QuantityExtension")
elseif (EQUATIONS STREQUAL "viscoelastic2")
  set(PARAMETER_MATERIAL "Viscoelastic")
  set(PARAMETER_VISCOMODE "AnelasticTensor")
else()
  capitalize(${EQUATIONS} PARAMETER_MATERIAL)
  set(PARAMETER_VISCOMODE "None")
endif()

capitalize(${DR_QUAD_RULE} PARAMETER_DRQUADRULE)

if (PRECISION STREQUAL "single")
  set(PARAMETER_REALTYPE "F32")
else()
  set(PARAMETER_REALTYPE "F64")
endif()

configure_file("Config.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/Config.h")

# Generate version.h
include(GetGitRevisionDescription)

# get GIT info
git_describe(PACKAGE_GIT_VERSION --tags --always --dirty=\ \(dirty\) --broken=\ \(broken\))
if (${PACKAGE_GIT_VERSION} MATCHES "NOTFOUND")
  set(PACKAGE_GIT_VERSION "(unknown)")
  set(PACKAGE_GIT_HASH "(unknown)")
  set(PACKAGE_GIT_TIMESTAMP "9999-12-31T00:00:00+00:00")
else()
  get_git_commit_info(PACKAGE_GIT_HASH PACKAGE_GIT_TIMESTAMP)
endif()
string(SUBSTRING ${PACKAGE_GIT_TIMESTAMP} 0 4 PACKAGE_GIT_YEAR)

# write file and print info
configure_file("Version.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/Version.h")
message(STATUS "Version: " ${PACKAGE_GIT_VERSION})
message(STATUS "Last commit: ${PACKAGE_GIT_HASH} at ${PACKAGE_GIT_TIMESTAMP}")
