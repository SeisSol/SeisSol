# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-kernels OBJECT
    DynamicRupture.cpp
    Plasticity.cpp
    PointSourceClusterOnDevice.cpp
    PointSourceClusterOnHost.cpp
    Receiver.cpp
    TimeCommon.cpp
    Touch.cpp

    LinearCK/GravitationalFreeSurfaceBC.cpp
)

target_link_libraries(seissol-kernels PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-kernels)
