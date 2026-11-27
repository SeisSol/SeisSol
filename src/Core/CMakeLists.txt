# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_subdirectory(Common)

add_subdirectory(Geometry)

add_subdirectory(IO)

add_subdirectory(Memory)

add_subdirectory(Modules)

add_subdirectory(Monitoring)

add_subdirectory(Numerical)

add_subdirectory(Parallel)

add_subdirectory(TimeStepping)

# The layers of the core (see scripts/ci/layering.py): a target uses only the ones below it; the headers
# of Core/Data, Core/Kernels and Core/Run come with the targets that include them.
target_link_libraries(seissol-parallel PUBLIC seissol-common)
target_link_libraries(seissol-geometry PUBLIC seissol-common seissol-parallel)
target_link_libraries(seissol-numerical PUBLIC seissol-common seissol-geometry seissol-parallel)
target_link_libraries(seissol-memory PUBLIC seissol-common seissol-parallel)
target_link_libraries(seissol-modules PUBLIC seissol-memory seissol-parallel)
target_link_libraries(seissol-io PUBLIC
    seissol-common seissol-memory seissol-modules seissol-numerical seissol-parallel)
target_link_libraries(seissol-monitoring PUBLIC
    seissol-common seissol-config seissol-io seissol-memory seissol-modules seissol-numerical
    seissol-parallel)
target_link_libraries(seissol-timestepping PUBLIC
    seissol-common seissol-geometry seissol-io seissol-memory seissol-monitoring
    seissol-numerical seissol-parallel)

# The targets of the core, with the description of the build (seissol-config). They link on their
# own, with the generated code: see seissol-core-link with the tests.
set_property(GLOBAL PROPERTY SEISSOL_CORE_TARGETS
    seissol-common
    seissol-config
    seissol-geometry
    seissol-io
    seissol-memory
    seissol-modules
    seissol-monitoring
    seissol-numerical
    seissol-parallel
    seissol-timestepping
)
