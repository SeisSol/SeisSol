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
    Datafield/Interpolation.cpp

    Scripting/CompiledReader.cpp
    Scripting/DataReader.cpp
    Scripting/DataTable.cpp
    Scripting/EasiReader.cpp
    Scripting/EasiWalker.cpp
    Scripting/LuaReader.cpp
    Scripting/LuaTracer.cpp
    Scripting/ReaderBuilder.cpp
)

# The samples of a grid are the reference a compiled kernel that samples it itself (Expr/RtcGpu.cpp)
# is checked against bit for bit: the same operations in the same order, without contraction of
# floating-point expressions, as the interpreter of src/Expr.
if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # -fp-model=fast, icpx's default from -O2 on, is -ffast-math
    set_source_files_properties(Datafield/Interpolation.cpp PROPERTIES COMPILE_OPTIONS
        "-ffp-contract=off;-fp-model=precise")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set_source_files_properties(Datafield/Interpolation.cpp PROPERTIES COMPILE_OPTIONS
        "-ffp-contract=off")
endif()

# seissol-external: lmathx, which the Lua reader loads as easi does (SEISSOL_LUA_LMATHX)
target_link_libraries(seissol-reader PUBLIC seissol-common-properties seissol-external)
target_link_libraries(seissol-lib PUBLIC seissol-reader)
