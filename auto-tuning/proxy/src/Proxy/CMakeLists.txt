# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-proxy-lib
    Allocator.cpp
    Common.cpp
    Kernel.cpp
    KernelHost.cpp
    KernelDevice.cpp
    Runner.cpp
    Tools.cpp
)

target_link_libraries(seissol-proxy-lib PUBLIC seissol-lib)
