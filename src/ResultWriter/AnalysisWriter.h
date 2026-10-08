// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_ANALYSISWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_ANALYSISWRITER_H_

#include "Geometry/MeshReader.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "Parallel/MPI.h"
#include "Physics/InitialField.h"

#include <array>
#include <cmath>

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::writer {
class AnalysisWriter {
  private:
  seissol::SeisSol& seissolInstance_;

  struct Data {
    double val;
    int rank;
  };

  bool isEnabled_{false}; // TODO(Lukas) Do we need this?
  const seissol::geometry::MeshReader* meshReader_{};

  std::string fileNamePrefix_;

  /// Compares the cells of the configuration `Cfg` with the initial condition; `configLabel` tells
  /// the configuration apart in the log, `fileName` is the table the errors go to.
  template <typename Cfg>
  void printAnalysisOf(double simulationTime,
                       const std::string& configLabel,
                       const std::string& fileName);

  public:
  explicit AnalysisWriter(seissol::SeisSol& seissolInstance) : seissolInstance_(seissolInstance) {}

  void init(const seissol::geometry::MeshReader* meshReader, std::string_view fileNamePrefix) {
    isEnabled_ = true;
    this->meshReader_ = meshReader;
    fileNamePrefix_ = std::string(fileNamePrefix);
  }

  void printAnalysis(double simulationTime);
}; // class AnalysisWriter

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_ANALYSISWRITER_H_
