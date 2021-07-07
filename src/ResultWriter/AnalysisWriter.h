#ifndef ANALYSISWRITER_H
#define ANALYSISWRITER_H

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

#include "Solver/Interoperability.h"
#include "Physics/InitialField.h"
#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/Transformation.h"
#include "Numerical_aux/BasisFunction.h"
#include "Parallel/MPI.h"
#include "Initializer/tree/Lut.hpp"

#include <Geometry/MeshReader.h>

extern seissol::Interoperability e_interoperability;

namespace seissol::writer {
class CsvAnalysisWriter {
public:
  CsvAnalysisWriter(std::string fileName);

  void writeHeader();

  void addObservation(std::string_view variable,
                      std::string_view normType,
                      real error);

  void enable();

  ~CsvAnalysisWriter();
private:
  std::ofstream out;
  bool isEnabled;
  std::string fileName;
};

  class AnalysisWriter {
private:
    struct data {
      double val;
      int rank;
    };

    bool isEnabled; // TODO(Lukas) Do we need this?
    const MeshReader* meshReader;

    std::string fileName;
public:
  AnalysisWriter() :
    isEnabled(false) { }

    void init(const MeshReader* meshReader,
              std::string_view fileNamePrefix) {
      isEnabled = true;
      this->meshReader = meshReader;
      fileName = std::string(fileNamePrefix) + "_analysis.csv";
    }  

    
    void printAnalysis(double simulationTime);
  }; // class AnalysisWriter


} // namespace seissol::writer
#endif // ANALYSISWRITER_H
