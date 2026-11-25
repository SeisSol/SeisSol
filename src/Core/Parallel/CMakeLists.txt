# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-parallel OBJECT
    Helper.cpp
    HelperThread.cpp
    MPI.cpp
    OpenMP.cpp
    ParallelFS.cpp
    Pin.cpp
    SystemInfo.cpp
    Trace.cpp

    Runtime/EventPool.cpp
    Runtime/HostActor.cpp
    Runtime/HostBufferPool.cpp
    Runtime/Stream.cpp
)

target_link_libraries(seissol-parallel PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-parallel)
