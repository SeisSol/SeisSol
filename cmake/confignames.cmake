# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The names of configurations, as seissol::configName writes them:
#   <material>[-m<mechanisms>]-o<order>-<f32|f64>[-s<fused simulations>]
# e.g. elastic-o4-f32 or viscoelastic-m3-o6-f64; and named sets of them, for CONFIGS and
# EXTRA_CONFIGS. Every part after the material carries a letter that says what it is. The solver of a
# material and the quadrature rule of the dynamic rupture faces are the same for all configurations
# of an executable (SOLVER and DR_QUAD_RULE), so they are not part of a name.

# The materials a configuration can have, by the name they have in a configuration name.
set(SEISSOL_MATERIALS elastic anisotropic viscoelastic poroelastic acoustic viscoacoustic)

# The named sets: SEISSOL_CONFIG_SETS lists their names, SEISSOL_CONFIG_SET_<name> the
# configurations of each.
#
# ci-cpu: what the CPU CI builds, on Linux and on macOS, in both precisions: the six materials in
# order 6, viscoelastic and viscoacoustic with 3 mechanisms, and elastic with 8 fused simulations in
# order 3, the order of the fused cases of precomputed-seissol. Building it into one executable takes
# a host architecture whose vectors the 8 fused simulations fill in both precisions, e.g. hsw or neon.
#
# ci-gpu: what the GPU CI builds with TensorForge, in both precisions: the six materials in order 6,
# viscoelastic and viscoacoustic with 3 mechanisms, and elastic with 32 fused simulations in order 3.
#
# ci-sanitize: what the sanitizer CI builds: the configurations of ci-cpu in double precision; in
# single precision, the same code runs with another type.
#
# ci-interaction: configurations whose cells are face neighbors across orders and precisions, for
# the conversions of their time integrals (codegen/kernels/configboundary.py): the solids lift the
# trace of the order 6 of their family to order 5 and pad the order 4 of the fluids, the fluids lift
# the order 6 of the solids to their order 4, and each family has a configuration in single and one
# in double precision; an anisotropic and an elastic solid of the same order meet as well, which
# pose the Riemann problem in different materials, and the anisotropic solid meets a fluid. Its unit
# tests, its verification cases with ConfigMap (interaction-* in tests/verification/Cases.cmake) and
# its sanitizer job check these conversions, which no other set takes, since all of their
# configurations of one simulation share their order.
#
# o2 to o8: one set per order, the six materials in that order and both precisions, viscoelastic and
# viscoacoustic with 3 mechanisms.
#
# all: the sets o2 to o8 together, o6 first: a run that names no configuration takes elastic-o6-f64,
# the configuration of a build with the default parameters.
set(SEISSOL_CONFIG_SETS ci-cpu ci-gpu ci-sanitize ci-interaction)

set(SEISSOL_CONFIG_SET_ci-cpu)
foreach(_precision f64 f32)
  list(APPEND SEISSOL_CONFIG_SET_ci-cpu
    elastic-o6-${_precision}
    acoustic-o6-${_precision}
    anisotropic-o6-${_precision}
    poroelastic-o6-${_precision}
    viscoelastic-m3-o6-${_precision}
    viscoacoustic-m3-o6-${_precision}
    elastic-o3-${_precision}-s8)
endforeach()

set(SEISSOL_CONFIG_SET_ci-gpu)
foreach(_precision f64 f32)
  list(APPEND SEISSOL_CONFIG_SET_ci-gpu
    elastic-o6-${_precision}
    acoustic-o6-${_precision}
    anisotropic-o6-${_precision}
    poroelastic-o6-${_precision}
    viscoelastic-m3-o6-${_precision}
    viscoacoustic-m3-o6-${_precision}
    elastic-o3-${_precision}-s32)
endforeach()

set(SEISSOL_CONFIG_SET_ci-sanitize)
foreach(_config IN LISTS SEISSOL_CONFIG_SET_ci-cpu)
  if (_config MATCHES "-f64")
    list(APPEND SEISSOL_CONFIG_SET_ci-sanitize ${_config})
  endif()
endforeach()

set(SEISSOL_CONFIG_SET_ci-interaction
  elastic-o6-f64
  anisotropic-o6-f64
  viscoelastic-m3-o5-f32
  acoustic-o4-f64
  viscoacoustic-m3-o4-f32)

foreach(_order 2 3 4 5 6 7 8)
  list(APPEND SEISSOL_CONFIG_SETS o${_order})
  set(SEISSOL_CONFIG_SET_o${_order})
  foreach(_precision f64 f32)
    list(APPEND SEISSOL_CONFIG_SET_o${_order}
      elastic-o${_order}-${_precision}
      acoustic-o${_order}-${_precision}
      anisotropic-o${_order}-${_precision}
      poroelastic-o${_order}-${_precision}
      viscoelastic-m3-o${_order}-${_precision}
      viscoacoustic-m3-o${_order}-${_precision})
  endforeach()
endforeach()

list(APPEND SEISSOL_CONFIG_SETS all)
set(SEISSOL_CONFIG_SET_all ${SEISSOL_CONFIG_SET_o6})
foreach(_order 2 3 4 5 7 8)
  list(APPEND SEISSOL_CONFIG_SET_all ${SEISSOL_CONFIG_SET_o${_order}})
endforeach()

# The name of a configuration.
function(seissol_config_name output material mechanisms order precision simulations)
  set(_name "${material}")
  if (mechanisms GREATER 0)
    string(APPEND _name "-m${mechanisms}")
  endif()
  if ("${precision}" STREQUAL "single")
    string(APPEND _name "-o${order}-f32")
  else()
    string(APPEND _name "-o${order}-f64")
  endif()
  if (NOT simulations EQUAL 1)
    string(APPEND _name "-s${simulations}")
  endif()
  set(${output} "${_name}" PARENT_SCOPE)
endfunction()

# Parses the configuration name `name` into <prefix>_MATERIAL, _MECHANISMS, _ORDER, _PRECISION
# (single or double) and _SIMULATIONS. Stops with an error that names `context` if `name` is no valid
# configuration, or not spelled as seissol::configName spells it. Needs the options of
# process_users_input.cmake.
function(seissol_parse_config_name prefix name context)
  # a spelling with a solver and a quadrature rule, or with mechanisms without their m: say which
  # name stands for the configuration
  if (name MATCHES
      "^([a-z]+)-(linearck|linearckanelastic|stp)(-m[0-9]+)?(-o[0-9]+-f(32|64))-[a-z]+(-s[0-9]+)?$")
    message(FATAL_ERROR "${context}: write ${name} as "
      "${CMAKE_MATCH_1}${CMAKE_MATCH_3}${CMAKE_MATCH_4}${CMAKE_MATCH_6}; the solver and the "
      "quadrature rule of the dynamic rupture faces are those of the build (SOLVER, DR_QUAD_RULE).")
  endif()
  if (name MATCHES "^([a-z]+)-([0-9]+)(-o[0-9]+-f(32|64)(-s[0-9]+)?)$")
    message(FATAL_ERROR "${context}: write ${name} as "
      "${CMAKE_MATCH_1}-m${CMAKE_MATCH_2}${CMAKE_MATCH_3}.")
  endif()
  if (NOT name MATCHES "^([a-z]+)(-m([0-9]+))?-o([0-9]+)-(f32|f64)(-s([0-9]+))?$")
    string(REPLACE ";" ", " _sets "${SEISSOL_CONFIG_SETS}")
    message(FATAL_ERROR "${context}: \"${name}\" is neither the name of a configuration, "
      "<material>[-m<mechanisms>]-o<order>-<f32|f64>[-s<fused simulations>], "
      "nor of a set of them (${_sets}).")
  endif()
  set(_material "${CMAKE_MATCH_1}")
  set(_mechanisms "${CMAKE_MATCH_3}")
  set(_order "${CMAKE_MATCH_4}")
  set(_precision "${CMAKE_MATCH_5}")
  set(_simulations "${CMAKE_MATCH_7}")
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

  check_parameter("The material of ${name}" ${_material} "${SEISSOL_MATERIALS}")
  check_parameter("The order of ${name}" ${_order} "${ORDER_OPTIONS}")
  if ((_material MATCHES "visco.?") AND (_mechanisms LESS 1))
    message(FATAL_ERROR "${context}: ${name} needs relaxation mechanisms (-m<number>).")
  endif()
  if ((NOT _material MATCHES "visco.?") AND (_mechanisms GREATER 0))
    message(FATAL_ERROR "${context}: ${_material} in ${name} has no relaxation mechanisms.")
  endif()

  seissol_config_name(_canonical ${_material} ${_mechanisms} ${_order} ${_precision} ${_simulations})
  if (NOT "${_canonical}" STREQUAL "${name}")
    message(FATAL_ERROR "${context}: write ${name} as ${_canonical}.")
  endif()

  set(${prefix}_MATERIAL ${_material} PARENT_SCOPE)
  set(${prefix}_MECHANISMS ${_mechanisms} PARENT_SCOPE)
  set(${prefix}_ORDER ${_order} PARENT_SCOPE)
  set(${prefix}_PRECISION ${_precision} PARENT_SCOPE)
  set(${prefix}_SIMULATIONS ${_simulations} PARENT_SCOPE)
endfunction()

# The C++ type of the configuration named `name`: its parts capitalized and joined, after Config,
# e.g. ConfigViscoelasticM3O6F64 for viscoelastic-m3-o6-f64.
function(seissol_config_type output name)
  set(_type Config)
  string(REPLACE "-" ";" _parts "${name}")
  foreach(_part IN LISTS _parts)
    string(SUBSTRING "${_part}" 0 1 _first)
    string(SUBSTRING "${_part}" 1 -1 _rest)
    string(TOUPPER "${_first}" _first)
    string(APPEND _type "${_first}${_rest}")
  endforeach()
  set(${output} ${_type} PARENT_SCOPE)
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

# The solver that advances the cells of `material` in this build: SOLVER, where the material can be
# advanced with it, and the first solver of the material otherwise (see process_users_input.cmake).
function(seissol_solver_of output material)
  set(_solvers ${SOLVERS_${material}})
  list(GET _solvers 0 _solver)
  if ((NOT SOLVER STREQUAL "auto") AND (SOLVER IN_LIST _solvers))
    set(_solver ${SOLVER})
  endif()
  set(${output} ${_solver} PARENT_SCOPE)
endfunction()
