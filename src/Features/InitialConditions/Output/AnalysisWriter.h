// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_INITIALCONDITIONS_OUTPUT_ANALYSISWRITER_H_
#define SEISSOL_SRC_FEATURES_INITIALCONDITIONS_OUTPUT_ANALYSISWRITER_H_

#include "Core/Geometry/MeshReader.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Numerical/BasisFunction.h"
#include "Core/Numerical/Quadrature.h"
#include "Core/Numerical/Transformation.h"
#include "Core/Parallel/MPI.h"
#include "Core/Run/Context.h"
#include "Features/InitialConditions/InitializationParameters.h"
#include "WavePropagation/InitialField.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace seissol::writer {
class AnalysisWriter {
  private:
  RunContext& run_;

  bool isEnabled_{false}; // TODO(Lukas) Do we need this?
  const seissol::geometry::MeshReader* meshReader_{};

  std::string fileNamePrefix_;

  // what init() gives: the initial conditions, as set up and per configuration, and the cells
  initializer::parameters::InitializationParameters initialization_;
  const std::vector<std::vector<std::unique_ptr<physics::InitialField>>>* initialConditions_{
      nullptr};
  LTS::Storage* storage_{nullptr};
  std::size_t simulationCount_{1};

  /// Compares the cells of the configuration `Cfg` with the initial condition; `configLabel` tells
  /// the configuration apart in the log, `fileName` is the table the errors go to.
  template <typename Cfg>
  void printAnalysisOf(double simulationTime,
                       const std::string& configLabel,
                       const std::string& fileName);

  public:
  explicit AnalysisWriter(RunContext& run) : run_(run) {}

  /// The cells of `storage` are compared with their `initialConditions` (by configuration), for
  /// `simulationCount` simulations.
  void init(
      const seissol::geometry::MeshReader* meshReader,
      std::string_view fileNamePrefix,
      const initializer::parameters::InitializationParameters& initialization,
      const std::vector<std::vector<std::unique_ptr<physics::InitialField>>>& initialConditions,
      LTS::Storage& storage,
      std::size_t simulationCount) {
    isEnabled_ = true;
    this->meshReader_ = meshReader;
    fileNamePrefix_ = std::string(fileNamePrefix);
    initialization_ = initialization;
    initialConditions_ = &initialConditions;
    storage_ = &storage;
    simulationCount_ = simulationCount;
  }

  void printAnalysis(double simulationTime);
}; // class AnalysisWriter

} // namespace seissol::writer

#endif // SEISSOL_SRC_FEATURES_INITIALCONDITIONS_OUTPUT_ANALYSISWRITER_H_
