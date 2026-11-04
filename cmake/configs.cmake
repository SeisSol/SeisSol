# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The configurations built into the executable, in the order of their ids. The first one is the one
# that EQUATIONS, NUMBER_OF_MECHANISMS, ORDER, PRECISION and NUMBER_OF_FUSED_SIMULATIONS describe;
# EXTRA_CONFIGS adds further ones by their name or by a named set of them (see
# cmake/confignames.cmake), e.g. elastic-o4-f32. With CONFIGS, its first configuration sets these
# variables (see cmake/process_users_input.cmake), and its others are the further ones.
#
# Sets, as lists with one entry per configuration:
#   SEISSOL_CONFIG_TYPES (the C++ type: Config0, Config1, ...), SEISSOL_CONFIG_MATERIALS,
#   SEISSOL_CONFIG_SOLVERS (the solver of the material in this build, see SOLVER),
#   SEISSOL_CONFIG_MECHANISMS, SEISSOL_CONFIG_ORDERS, SEISSOL_CONFIG_PRECISIONS (single or double),
#   SEISSOL_CONFIG_SIMULATIONS, SEISSOL_CONFIG_NAMES (e.g. elastic-o6-f64);
# and SEISSOL_SOLVERS_USED, the solvers of all configurations, each once.

set(EXTRA_CONFIGS "" CACHE STRING
  "Further configurations to build into the executable, by name (e.g. elastic-o4-f32)")

seissol_expand_configs(_further_configs "${EXTRA_CONFIGS}")
set(_further_context "EXTRA_CONFIGS")
if (NOT "${CONFIGS}" STREQUAL "")
  if (NOT "${EXTRA_CONFIGS}" STREQUAL "")
    message(FATAL_ERROR "Give the configurations either as CONFIGS or as EXTRA_CONFIGS (with "
      "EQUATIONS, ORDER, ...), not both.")
  endif()
  set(_further_configs ${SEISSOL_CONFIGS_FURTHER})
  set(_further_context "CONFIGS")
endif()

set(SEISSOL_CONFIG_TYPES Config0)
set(SEISSOL_CONFIG_MATERIALS ${EQUATIONS})
set(SEISSOL_CONFIG_MECHANISMS ${NUMBER_OF_MECHANISMS})
set(SEISSOL_CONFIG_ORDERS ${ORDER})
set(SEISSOL_CONFIG_PRECISIONS ${PRECISION})
set(SEISSOL_CONFIG_SIMULATIONS ${NUMBER_OF_FUSED_SIMULATIONS})
seissol_config_name(_config_names ${EQUATIONS} ${NUMBER_OF_MECHANISMS} ${ORDER} ${PRECISION}
  ${NUMBER_OF_FUSED_SIMULATIONS})

set(_config_index 0)
foreach(_name IN LISTS _further_configs)
  seissol_parse_config_name(_config ${_name} ${_further_context})
  if (_name IN_LIST _config_names)
    message(FATAL_ERROR "${_further_context}: ${_name} is built already.")
  endif()
  list(APPEND _config_names ${_name})

  math(EXPR _config_index "${_config_index} + 1")
  list(APPEND SEISSOL_CONFIG_TYPES Config${_config_index})
  list(APPEND SEISSOL_CONFIG_MATERIALS ${_config_MATERIAL})
  list(APPEND SEISSOL_CONFIG_MECHANISMS ${_config_MECHANISMS})
  list(APPEND SEISSOL_CONFIG_ORDERS ${_config_ORDER})
  list(APPEND SEISSOL_CONFIG_PRECISIONS ${_config_PRECISION})
  list(APPEND SEISSOL_CONFIG_SIMULATIONS ${_config_SIMULATIONS})
endforeach()

list(LENGTH SEISSOL_CONFIG_TYPES SEISSOL_CONFIG_COUNT)
math(EXPR SEISSOL_CONFIG_LAST "${SEISSOL_CONFIG_COUNT} - 1")

# the solver of each configuration follows from its material, as SOLVER chooses it for the build
set(SEISSOL_CONFIG_SOLVERS)
foreach(_material IN LISTS SEISSOL_CONFIG_MATERIALS)
  seissol_solver_of(_solver ${_material})
  list(APPEND SEISSOL_CONFIG_SOLVERS ${_solver})
endforeach()

set(SEISSOL_SOLVERS_USED ${SEISSOL_CONFIG_SOLVERS})
list(REMOVE_DUPLICATES SEISSOL_SOLVERS_USED)

if ((NOT SOLVER STREQUAL "auto") AND (NOT SOLVER IN_LIST SEISSOL_SOLVERS_USED))
  set(_materials ${SEISSOL_CONFIG_MATERIALS})
  list(REMOVE_DUPLICATES _materials)
  if (SOLVER STREQUAL "stp" AND "${_materials}" MATCHES "visco")
    message(FATAL_ERROR
      "SOLVER=stp cannot be combined with ${_materials} yet: the space-time predictor substitutes "
      "one entry per stiff source row, and relaxation needs several per row.")
  endif()
  message(FATAL_ERROR "SOLVER=${SOLVER}: no material of the build (${_materials}) can be advanced "
    "with it.")
endif()

set(SEISSOL_CONFIG_NAMES ${_config_names})

message(STATUS "Configurations: ${_config_names}")
set(_material_solvers)
foreach(_material _solver IN ZIP_LISTS SEISSOL_CONFIG_MATERIALS SEISSOL_CONFIG_SOLVERS)
  list(APPEND _material_solvers "${_material}: ${_solver}")
endforeach()
list(REMOVE_DUPLICATES _material_solvers)
string(REPLACE ";" ", " _material_solvers "${_material_solvers}")
message(STATUS "Solvers: ${_material_solvers}")
message(STATUS "Dynamic rupture quadrature rule: ${DR_QUAD_RULE}")
