// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_RESULTWRITER_ANALYSISWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_ANALYSISWRITER_H_

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

#include "Initializer/Tree/Lut.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "Parallel/MPI.h"
#include "Physics/InitialField.h"

#include "Geometry/MeshReader.h"

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::writer {
class CsvAnalysisWriter {
  public:
  CsvAnalysisWriter(std::string fileName);

  void writeHeader();

  void addObservation(std::string_view variable, std::string_view normType, real error);

  void enable();

  ~CsvAnalysisWriter();

  private:
  std::ofstream out;
  bool isEnabled;
  std::string fileName;
};

class AnalysisWriter {
  private:
  seissol::SeisSol& seissolInstance;

  struct Data {
    double val;
    int rank;
  };

  bool isEnabled; // TODO(Lukas) Do we need this?
  const seissol::geometry::MeshReader* meshReader;

  std::string fileName;

  public:
  AnalysisWriter(seissol::SeisSol& seissolInstance)
      : seissolInstance(seissolInstance), isEnabled(false) {}

  void init(const seissol::geometry::MeshReader* meshReader, std::string_view fileNamePrefix) {
    isEnabled = true;
    this->meshReader = meshReader;
    fileName = std::string(fileNamePrefix) + "-analysis.csv";
  }

  void printAnalysis(double simulationTime);
}; // class AnalysisWriter

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_ANALYSISWRITER_H_
