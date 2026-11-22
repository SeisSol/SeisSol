# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-dynamicrupture OBJECT
    Misc.cpp
    Output/Builders/ReceiverBasedOutputBuilder.cpp
    Output/FaultRefiner/FaultRefiners.cpp
    Output/OutputAux.cpp
    Output/OutputManager.cpp
    Output/PickpointOutput.cpp
    Output/ReceiverBasedOutput.cpp
    Setup/DynamicRuptureMatrices.cpp
    Setup/FaultParameterDB.cpp
    Setup/Factory.cpp
    Setup/Initializer/BaseDRInitializer.cpp
    Setup/Initializer/ImposedSlipRatesInitializer.cpp
    Setup/Initializer/LinearSlipWeakeningInitializer.cpp
    Setup/Initializer/RateAndStateInitializer.cpp
    Solver/DynamicRupture.cpp
    Solver/DynamicRuptureCluster.cpp
    Solver/DynamicRuptureClusterFactory.cpp
    Solver/DynamicRuptureFeature.cpp
    Solver/DynamicRuptureRecorder.cpp
    Solver/FrictionLaws/CpuImpl/SourceTimeFunction.cpp
    Solver/FrictionLaws/CpuImpl/ThermalPressurization/ThermalPressurization.cpp
    Solver/FrictionLaws/FrictionSolver.cpp
)

target_link_libraries(seissol-dynamicrupture PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-dynamicrupture)
