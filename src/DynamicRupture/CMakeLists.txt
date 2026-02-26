# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-dynamicrupture OBJECT
    Factory.cpp
    Misc.cpp
    FrictionLaws/FrictionSolver.cpp
    FrictionLaws/CpuImpl/LinearSlipWeakening.cpp
    FrictionLaws/CpuImpl/NoFault.cpp
    FrictionLaws/CpuImpl/SourceTimeFunction.cpp
    FrictionLaws/CpuImpl/ThermalPressurization/ThermalPressurization.cpp
    FrictionLaws/GpuImpl/FrictionSolverDetails.cpp
    Initializer/BaseDRInitializer.cpp
    Initializer/ImposedSlipRatesInitializer.cpp
    Initializer/LinearSlipWeakeningInitializer.cpp
    Initializer/RateAndStateInitializer.cpp
    Output/OutputAux.cpp
    Output/OutputManager.cpp
    Output/ReceiverBasedOutput.cpp
    Output/Builders/ReceiverBasedOutputBuilder.cpp
    Output/FaultRefiner/FaultRefiners.cpp
)

target_link_libraries(seissol-dynamicrupture PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-dynamicrupture)
