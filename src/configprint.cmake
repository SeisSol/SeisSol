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

# the configurations built into the executable (cmake/configs.cmake): a type each, and the lists of
# them for Config.h
set(SEISSOL_CONFIG_STRUCTS "")
set(SEISSOL_CONFIG_LIST "")
set(SEISSOL_CONFIG_TYPE_NAMES "")
# the configurations each solver advances, for the instantiation of its kernels
set(SEISSOL_CONFIGS_LINEARCK "")
set(SEISSOL_CONFIGS_LINEARCKANELASTIC "")
set(SEISSOL_CONFIGS_STP "")
foreach(IDX RANGE ${SEISSOL_CONFIG_LAST})
  foreach(FIELD TYPES MATERIALS SOLVERS MECHANISMS ORDERS PRECISIONS DRQUADRULES SIMULATIONS)
    list(GET SEISSOL_CONFIG_${FIELD} ${IDX} CONFIG_${FIELD})
  endforeach()

  capitalize(${CONFIG_MATERIALS} PARAMETER_MATERIAL)
  capitalize(${CONFIG_DRQUADRULES} PARAMETER_DRQUADRULE)
  if (CONFIG_SOLVERS STREQUAL "linearck")
    set(PARAMETER_SOLVER "LinearCK")
  elseif (CONFIG_SOLVERS STREQUAL "linearckanelastic")
    set(PARAMETER_SOLVER "LinearCKAnelastic")
  elseif (CONFIG_SOLVERS STREQUAL "stp")
    set(PARAMETER_SOLVER "STP")
  else()
    message(FATAL_ERROR "Invalid solver: ${CONFIG_SOLVERS}")
  endif()
  if (CONFIG_PRECISIONS STREQUAL "single")
    set(PARAMETER_REALTYPE "F32")
  else()
    set(PARAMETER_REALTYPE "F64")
  endif()

  string(APPEND SEISSOL_CONFIG_STRUCTS
    "struct ${CONFIG_TYPES}\n"
    "    : ConfigOf<${CONFIG_ORDERS},\n"
    "               ${CONFIG_MECHANISMS},\n"
    "               model::MaterialType::${PARAMETER_MATERIAL},\n"
    "               RealType::${PARAMETER_REALTYPE},\n"
    "               SolverType::${PARAMETER_SOLVER},\n"
    "               DRQuadRuleType::${PARAMETER_DRQUADRULE},\n"
    "               ${CONFIG_SIMULATIONS}> {};\n")
  string(APPEND SEISSOL_CONFIG_LIST " X(::seissol::${CONFIG_TYPES})")
  list(APPEND SEISSOL_CONFIG_TYPE_NAMES "::seissol::${CONFIG_TYPES}")
  string(TOUPPER ${PARAMETER_SOLVER} PARAMETER_SOLVER_UPPER)
  string(APPEND SEISSOL_CONFIGS_${PARAMETER_SOLVER_UPPER} " X(::seissol::${CONFIG_TYPES})")
endforeach()
string(JOIN ", " SEISSOL_CONFIG_TYPE_LIST ${SEISSOL_CONFIG_TYPE_NAMES})

configure_file("Config.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/Config.h")

# Generate BuildInfo.cpp
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
configure_file("BuildInfo.cpp.in"
               "${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.cpp")
message(STATUS "Version: " ${PACKAGE_GIT_VERSION})
message(STATUS "Last commit: ${PACKAGE_GIT_HASH} at ${PACKAGE_GIT_TIMESTAMP}")

add_library(seissol-config OBJECT
${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.cpp
)

target_link_libraries(seissol-config PUBLIC seissol-common-properties)
