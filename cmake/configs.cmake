# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The configurations built into the executable, in the order of their ids. The first one is the one
# that EQUATIONS, SOLVER, NUMBER_OF_MECHANISMS, ORDER, PRECISION, DR_QUAD_RULE and
# NUMBER_OF_FUSED_SIMULATIONS describe; EXTRA_CONFIGS adds further ones by their name, as
# seissol::configName writes it: <material>-<solver>[-m<mechanisms>]-o<order>-<f32|f64>-<dr quadrature
# rule>[-s<fused simulations>], e.g. elastic-linearck-o4-f32-stroud.
#
# Sets, as lists with one entry per configuration:
#   SEISSOL_CONFIG_TYPES (the C++ type: Config, Config1, ...), SEISSOL_CONFIG_MATERIALS,
#   SEISSOL_CONFIG_SOLVERS, SEISSOL_CONFIG_MECHANISMS, SEISSOL_CONFIG_ORDERS,
#   SEISSOL_CONFIG_PRECISIONS (single or double), SEISSOL_CONFIG_DRQUADRULES,
#   SEISSOL_CONFIG_SIMULATIONS;
# and SEISSOL_SOLVERS_USED, the solvers of all configurations, each once.

set(EXTRA_CONFIGS "" CACHE STRING
  "Further configurations to build into the executable, by name (e.g. elastic-linearck-o4-f32-stroud)")

# The name of a configuration, as seissol::configName writes it.
function(seissol_config_name output material solver mechanisms order precision drquadrule simulations)
  set(_name "${material}-${solver}")
  if (mechanisms GREATER 0)
    string(APPEND _name "-m${mechanisms}")
  endif()
  if ("${precision}" STREQUAL "single")
    string(APPEND _name "-o${order}-f32-${drquadrule}")
  else()
    string(APPEND _name "-o${order}-f64-${drquadrule}")
  endif()
  if (NOT simulations EQUAL 1)
    string(APPEND _name "-s${simulations}")
  endif()
  set(${output} "${_name}" PARENT_SCOPE)
endfunction()

set(SEISSOL_CONFIG_TYPES Config)
set(SEISSOL_CONFIG_MATERIALS ${EQUATIONS})
set(SEISSOL_CONFIG_SOLVERS ${SOLVER})
set(SEISSOL_CONFIG_MECHANISMS ${NUMBER_OF_MECHANISMS})
set(SEISSOL_CONFIG_ORDERS ${ORDER})
set(SEISSOL_CONFIG_PRECISIONS ${PRECISION})
set(SEISSOL_CONFIG_DRQUADRULES ${DR_QUAD_RULE})
set(SEISSOL_CONFIG_SIMULATIONS ${NUMBER_OF_FUSED_SIMULATIONS})
seissol_config_name(_config_names ${EQUATIONS} ${SOLVER} ${NUMBER_OF_MECHANISMS} ${ORDER} ${PRECISION}
  ${DR_QUAD_RULE} ${NUMBER_OF_FUSED_SIMULATIONS})

set(_config_index 0)
foreach(_name IN LISTS EXTRA_CONFIGS)
  if (NOT _name MATCHES "^([a-z]+)-([a-z]+)(-m([0-9]+))?-o([0-9]+)-(f32|f64)-([a-z]+)(-s([0-9]+))?$")
    message(FATAL_ERROR "EXTRA_CONFIGS: \"${_name}\" is not the name of a configuration, "
      "<material>-<solver>[-m<mechanisms>]-o<order>-<f32|f64>-<dr quadrature rule>[-s<fused simulations>].")
  endif()
  set(_material "${CMAKE_MATCH_1}")
  set(_solver "${CMAKE_MATCH_2}")
  set(_mechanisms "${CMAKE_MATCH_4}")
  set(_order "${CMAKE_MATCH_5}")
  set(_precision "${CMAKE_MATCH_6}")
  set(_drquadrule "${CMAKE_MATCH_7}")
  set(_simulations "${CMAKE_MATCH_9}")
  if ("${_mechanisms}" STREQUAL "")
    set(_mechanisms 0)
  endif()
  if ("${_simulations}" STREQUAL "")
    set(_simulations 1)
  endif()
  if ("${_precision}" STREQUAL "f32")
    set(_precision single)
  else()
    set(_precision double)
  endif()

  check_parameter("The material of ${_name}" ${_material} "${EQUATIONS_OPTIONS}")
  if (NOT _solver IN_LIST SOLVERS_${_material})
    message(FATAL_ERROR "EXTRA_CONFIGS: ${_name} advances ${_material} with ${_solver}; "
      "available: ${SOLVERS_${_material}}.")
  endif()
  check_parameter("The order of ${_name}" ${_order} "${ORDER_OPTIONS}")
  check_parameter("The dynamic rupture quadrature rule of ${_name}" ${_drquadrule}
    "${DR_QUAD_RULE_OPTIONS}")
  if ((_material MATCHES "visco.?") AND (_mechanisms LESS 1))
    message(FATAL_ERROR "EXTRA_CONFIGS: ${_name} needs relaxation mechanisms (-m<number>).")
  endif()
  if ((NOT _material MATCHES "visco.?") AND (_mechanisms GREATER 0))
    message(FATAL_ERROR "EXTRA_CONFIGS: ${_material} in ${_name} has no relaxation mechanisms.")
  endif()

  seissol_config_name(_canonical ${_material} ${_solver} ${_mechanisms} ${_order} ${_precision}
    ${_drquadrule} ${_simulations})
  if (NOT "${_canonical}" STREQUAL "${_name}")
    message(FATAL_ERROR "EXTRA_CONFIGS: write ${_name} as ${_canonical}.")
  endif()
  if (_canonical IN_LIST _config_names)
    message(FATAL_ERROR "EXTRA_CONFIGS: ${_name} is built already.")
  endif()
  list(APPEND _config_names ${_canonical})

  math(EXPR _config_index "${_config_index} + 1")
  list(APPEND SEISSOL_CONFIG_TYPES Config${_config_index})
  list(APPEND SEISSOL_CONFIG_MATERIALS ${_material})
  list(APPEND SEISSOL_CONFIG_SOLVERS ${_solver})
  list(APPEND SEISSOL_CONFIG_MECHANISMS ${_mechanisms})
  list(APPEND SEISSOL_CONFIG_ORDERS ${_order})
  list(APPEND SEISSOL_CONFIG_PRECISIONS ${_precision})
  list(APPEND SEISSOL_CONFIG_DRQUADRULES ${_drquadrule})
  list(APPEND SEISSOL_CONFIG_SIMULATIONS ${_simulations})
endforeach()

list(LENGTH SEISSOL_CONFIG_TYPES SEISSOL_CONFIG_COUNT)
math(EXPR SEISSOL_CONFIG_LAST "${SEISSOL_CONFIG_COUNT} - 1")

set(SEISSOL_SOLVERS_USED ${SEISSOL_CONFIG_SOLVERS})
list(REMOVE_DUPLICATES SEISSOL_SOLVERS_USED)

message(STATUS "Configurations: ${_config_names}")
