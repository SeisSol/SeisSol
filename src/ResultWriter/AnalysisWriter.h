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

  std::string fileName_;

  public:
  explicit AnalysisWriter(seissol::SeisSol& seissolInstance) : seissolInstance_(seissolInstance) {}

  void init(const seissol::geometry::MeshReader* meshReader, std::string_view fileNamePrefix) {
    isEnabled_ = true;
    this->meshReader_ = meshReader;
    fileName_ = std::string(fileNamePrefix) + "-analysis.csv";
  }

  void printAnalysis(double simulationTime);
}; // class AnalysisWriter

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_ANALYSISWRITER_H_
