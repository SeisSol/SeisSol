#include "FaultRefiner.hpp"


namespace seissol::dr::output {
  std::unique_ptr<FaultRefinerInterface> getRefiner(const int strategy) {
    switch (strategy) {
      case 1:
        return std::unique_ptr<FaultRefinerInterface>(new TripleFaultFaceRefiner);
      case 2:
        return std::unique_ptr<FaultRefinerInterface>(new QuadFaultFaceRefiner);
      default:
        throw std::runtime_error("Unknown refinement strategy for Fault Face Refiner");
    }
  }
}