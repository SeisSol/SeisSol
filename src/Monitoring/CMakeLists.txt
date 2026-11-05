# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-monitoring OBJECT
    ActorStateStatistics.cpp
    Counters.cpp
    Energy/Discovery.cpp
    Energy/Library.cpp
    Energy/Sensor.cpp
    Energy/Sysfs.cpp
    Energy/Vendors.cpp
    EpochSynchronization.cpp
    FlopCounter.cpp
    LoopStatistics.cpp
    Memory.cpp
    RunSummary.cpp
    Stopwatch.cpp
    Unit.cpp
)

target_link_libraries(seissol-monitoring PUBLIC seissol-common-properties)
# the management libraries of the GPU vendors are opened at runtime (dlopen is in libdl before
# glibc 2.34)
target_link_libraries(seissol-monitoring PUBLIC ${CMAKE_DL_LIBS})
target_link_libraries(seissol-lib PUBLIC seissol-monitoring)
