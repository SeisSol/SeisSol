# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The names of configurations, as seissol::configName writes them:
#   <material>-<solver>[-m<mechanisms>]-o<order>-<f32|f64>-<dr quadrature rule>[-s<fused simulations>]
# e.g. elastic-linearck-o4-f32-stroud; and named sets of them, for CONFIGS and EXTRA_CONFIGS.

# The named sets: SEISSOL_CONFIG_SETS lists their names, SEISSOL_CONFIG_SET_<name> the
# configurations of each.
#
# ci-cpu: what the CPU CI builds, in order 6 and both precisions: the six materials, viscoelastic
# and viscoacoustic with their default solver and 3 mechanisms, and elastic with 8 fused
# simulations. Building it into one executable takes a host architecture whose vectors the 8 fused
# simulations fill in both precisions, e.g. hsw.
set(SEISSOL_CONFIG_SETS ci-cpu)

set(SEISSOL_CONFIG_SET_ci-cpu)
foreach(_precision f64 f32)
  list(APPEND SEISSOL_CONFIG_SET_ci-cpu
    elastic-linearck-o6-${_precision}-stroud
    acoustic-linearck-o6-${_precision}-stroud
    anisotropic-linearck-o6-${_precision}-stroud
    poroelastic-stp-o6-${_precision}-stroud
    viscoelastic-linearckanelastic-m3-o6-${_precision}-stroud
    viscoacoustic-linearckanelastic-m3-o6-${_precision}-stroud
    elastic-linearck-o6-${_precision}-stroud-s8)
endforeach()

# The name of a configuration.
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

# Parses the configuration name `name` into <prefix>_MATERIAL, _SOLVER, _MECHANISMS, _ORDER,
# _PRECISION (single or double), _DRQUADRULE and _SIMULATIONS. Stops with an error that names
# `context` if `name` is no valid configuration, or not spelled as seissol::configName spells it.
# Needs the options of process_users_input.cmake.
function(seissol_parse_config_name prefix name context)
  if (NOT name MATCHES "^([a-z]+)-([a-z]+)(-m([0-9]+))?-o([0-9]+)-(f32|f64)-([a-z]+)(-s([0-9]+))?$")
    message(FATAL_ERROR "${context}: \"${name}\" is neither the name of a configuration, "
      "<material>-<solver>[-m<mechanisms>]-o<order>-<f32|f64>-<dr quadrature rule>[-s<fused simulations>], "
      "nor of a set of them (${SEISSOL_CONFIG_SETS}).")
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

  check_parameter("The material of ${name}" ${_material} "${EQUATIONS_OPTIONS}")
  if (NOT _solver IN_LIST SOLVERS_${_material})
    message(FATAL_ERROR "${context}: ${name} advances ${_material} with ${_solver}; "
      "available: ${SOLVERS_${_material}}.")
  endif()
  check_parameter("The order of ${name}" ${_order} "${ORDER_OPTIONS}")
  check_parameter("The dynamic rupture quadrature rule of ${name}" ${_drquadrule}
    "${DR_QUAD_RULE_OPTIONS}")
  if ((_material MATCHES "visco.?") AND (_mechanisms LESS 1))
    message(FATAL_ERROR "${context}: ${name} needs relaxation mechanisms (-m<number>).")
  endif()
  if ((NOT _material MATCHES "visco.?") AND (_mechanisms GREATER 0))
    message(FATAL_ERROR "${context}: ${_material} in ${name} has no relaxation mechanisms.")
  endif()

  seissol_config_name(_canonical ${_material} ${_solver} ${_mechanisms} ${_order} ${_precision}
    ${_drquadrule} ${_simulations})
  if (NOT "${_canonical}" STREQUAL "${name}")
    message(FATAL_ERROR "${context}: write ${name} as ${_canonical}.")
  endif()

  set(${prefix}_MATERIAL ${_material} PARENT_SCOPE)
  set(${prefix}_SOLVER ${_solver} PARENT_SCOPE)
  set(${prefix}_MECHANISMS ${_mechanisms} PARENT_SCOPE)
  set(${prefix}_ORDER ${_order} PARENT_SCOPE)
  set(${prefix}_PRECISION ${_precision} PARENT_SCOPE)
  set(${prefix}_DRQUADRULE ${_drquadrule} PARENT_SCOPE)
  set(${prefix}_SIMULATIONS ${_simulations} PARENT_SCOPE)
endfunction()

# The configuration names that the list `names` stands for: its names, with each named set replaced
# by the configurations it holds.
function(seissol_expand_configs output names)
  set(_expanded)
  foreach(_name IN LISTS names)
    if (_name IN_LIST SEISSOL_CONFIG_SETS)
      list(APPEND _expanded ${SEISSOL_CONFIG_SET_${_name}})
    else()
      list(APPEND _expanded ${_name})
    endif()
  endforeach()
  set(${output} ${_expanded} PARENT_SCOPE)
endfunction()
