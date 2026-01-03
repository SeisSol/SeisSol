# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

add_library(seissol-resultwriter OBJECT
    AnalysisWriter.cpp
    AsyncIO.cpp
    ClusteringWriter.cpp
    EnergyOutput.cpp
    FaultWriter.cpp
    FaultWriterExecutor.cpp
    FreeSurfaceWriter.cpp
    FreeSurfaceWriterExecutor.cpp
    MiniSeisSolWriter.cpp
    PostProcessor.cpp
    ReceiverWriter.cpp
    ThreadsPinningWriter.cpp
    WaveFieldWriter.cpp
)

target_link_libraries(seissol-resultwriter PUBLIC seissol-common-properties)
target_link_libraries(seissol-lib PUBLIC seissol-resultwriter)
