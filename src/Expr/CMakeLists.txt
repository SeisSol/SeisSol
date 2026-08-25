# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-expr OBJECT
    Ir.cpp
    Program.cpp
    Binding.cpp
    Backend.cpp
    RtcCpu.cpp
    Validate.cpp
    Lower.cpp
    Interp.cpp
    SderivFrontend.cpp
)

target_link_libraries(seissol-expr PUBLIC seissol-common-properties)

# The compiled CPU backend loads its artifact with dlopen. It degrades to the
# interpreter at run time when the toolchain is unusable, so this is the only
# build-time requirement it has.
target_link_libraries(seissol-expr PUBLIC ${CMAKE_DL_LIBS})
target_link_libraries(seissol-lib PUBLIC seissol-expr)
