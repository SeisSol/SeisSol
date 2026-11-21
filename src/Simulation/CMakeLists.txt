# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-simulation OBJECT
    BuiltConfigs.cpp
    EnergyOutput.cpp
    Estimator.cpp
    InitProcedure/Init.cpp
    InitProcedure/InitIO.cpp
    InitProcedure/InitLayout.cpp
    InitProcedure/InitMesh.cpp
    InitProcedure/InitModel.cpp
    InitProcedure/InitSideConditions.cpp
    InitProcedure/Internal/Boundary.cpp
    InitProcedure/Internal/Buckets.cpp
    InitProcedure/Internal/ConfigBoundaryCheck.cpp
    InitProcedure/Internal/FaceTypeCheck.cpp
    InitProcedure/Internal/LtsSetup.cpp
    InitProcedure/Internal/MeshLayout.cpp
    InitProcedure/Internal/Recording.cpp
    InitProcedure/Internal/Scratchpads.cpp
    MemoryManager.cpp
    Parameters/CubeGeneratorParameters.cpp
    Parameters/DRParameters.cpp
    Parameters/InitializationParameters.cpp
    Parameters/LtsParameters.cpp
    Parameters/MeshParameters.cpp
    Parameters/ModelParameters.cpp
    Parameters/OutputParameters.cpp
    Parameters/ParameterOverrides.cpp
    Parameters/ParameterReader.cpp
    Parameters/SeisSolParameters.cpp
    Parameters/SourceParameters.cpp
    Preflight.cpp
    PreflightMesh.cpp
    Reader/AsagiModule.cpp
    Reader/AsagiReader.cpp
    Reader/EasiReader.cpp
    Reader/ReaderBuilder.cpp
    Simulator.cpp
)

target_link_libraries(seissol-simulation PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-simulation)
