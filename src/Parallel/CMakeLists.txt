# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-parallel OBJECT
    AcceleratorDevice.cpp
    DataCollector.cpp
    HelperThread.cpp
    MPI.cpp
    OpenMP.cpp
    Pin.cpp

    Runtime/Stream.cpp
)

target_link_libraries(seissol-parallel PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-parallel)
