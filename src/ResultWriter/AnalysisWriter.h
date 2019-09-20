#ifndef ANALYSISWRITER_H
#define ANALYSISWRITER_H

// TODO(Lukas) Clean up includes.
#include <array>
#include <cmath>

#include "Solver/Interoperability.h"
#include "Physics/InitialField.h"
#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/Transformation.h"
#include "Numerical_aux/BasisFunction.h"
#include "Parallel/MPI.h"
#include "Initializer/tree/Lut.hpp"

#include <Geometry/MeshReader.h>

//#include <Initializer/LTS.h>
//#include <Initializer/tree/LTSTree.hpp>

extern seissol::Interoperability e_interoperability;

namespace seissol {
namespace writer {
  class AnalysisWriter {
private:
    struct data {
      double val;
      int rank;
    };

    bool isEnabled; // TODO(Lukas) Do we need this?
    const MeshReader* meshReader;
public:
  AnalysisWriter() :
    isEnabled(false) { }

    void init(const MeshReader* meshReader) {
      isEnabled = true;
      this->meshReader = meshReader;
    }  

    
    void printAnalysis(double simulationTime);
  }; // class AnalysisWriter
} // namespace Writer
} // namespace Solver
#endif // ANALYSISWRITER_H
