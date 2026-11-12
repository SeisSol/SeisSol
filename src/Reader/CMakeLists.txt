# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-reader OBJECT
    Datafield/AsagiLite.cpp
    Datafield/AsagiModule.cpp
    Datafield/AsagiReader.cpp
    Datafield/Grid.cpp
    Datafield/GridFactory.cpp
    Datafield/GridUpdateModule.cpp
    Datafield/InlineGrid.cpp
    Datafield/ScecGrid.cpp
    Datafield/Interpolation.cpp

    Scripting/CompiledReader.cpp
    Scripting/DataReader.cpp
    Scripting/DataTable.cpp
    Scripting/EasiReader.cpp
    Scripting/EasiSpecial.cpp
    Scripting/EasiWalker.cpp
    Scripting/ImpalaToLua.cpp
    Scripting/LuaFrontend.cpp
    Scripting/LuaSyntax.cpp
    Scripting/ReaderBuilder.cpp
)

# Without contraction of floating-point expressions, and IEEE: what the reader computes is the bits
# of another implementation. The samples of a grid are the reference a compiled kernel that samples
# it itself (Expr/RtcGpu.cpp) is checked against bit for bit, as is the interpreter of src/Expr; the
# Lua frontend computes what the program does not as Lua does, and the easi walker what it computes
# itself (the coefficients of a layered model, the constants of a special component) as easi does.
# A contraction a*b + c would be a fused multiply-add here and not there.
if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # -fp-model=fast, icpx's default from -O2 on, is -ffast-math
    target_compile_options(seissol-reader PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off> $<$<COMPILE_LANGUAGE:CXX>:-fp-model=precise>)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(seissol-reader PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)
endif()

target_link_libraries(seissol-reader PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-reader)
