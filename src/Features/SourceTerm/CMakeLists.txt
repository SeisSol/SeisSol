# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-sourceterm OBJECT
    Physics/PointSource.cpp
    Setup/FSRMReader.cpp
    Setup/Manager.cpp
    Solver/PointSourceClusterOnHost.cpp
)

# the point sources on the device, whose kernels are in seissol-device-lib
if (WITH_GPU)
    target_sources(seissol-sourceterm PRIVATE Solver/PointSourceClusterOnDevice.cpp)
endif()

if (NETCDF)
    target_sources(seissol-sourceterm PRIVATE Setup/NRFReader.cpp)
endif()

target_link_libraries(seissol-sourceterm PUBLIC seissol-common-properties)
# a feature uses the wave propagation, and the core below it
target_link_libraries(seissol-sourceterm PUBLIC seissol-wavepropagation)
target_link_libraries(seissol-lib PUBLIC seissol-sourceterm)
