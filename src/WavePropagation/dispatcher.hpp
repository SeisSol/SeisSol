#pragma once

#include <memory>
#include "Initializer/LTS.h"
#include "Initializer/tree/Layer.hpp"

namespace seissol::waveprop {
class WavePropDispatcherBase {
  public:
  virtual void dispatchPredict(double timeStepSize, double correctionTime, bool resetBuffers) = 0;
  virtual void computePredictFlops(long long int& flopsNonZero, long long int& flopsHardware) = 0;
  virtual void dispatchNeighborCorrect(double timeStepSize, double subTimeStart) = 0;
  virtual void computeCorrectFlops(long long int& flopsNonZero,
                                   long long int& flopsHardware,
                                   long long int& drFlopsNonZero,
                                   long long int& drFlopsHardware) = 0;
  virtual void setTV(double tv) = 0;

  private:
};

std::unique_ptr<WavePropDispatcherBase> getDispatcher(const seissol::initializers::LTS& lts,
                                                      seissol::initializers::Layer& layer,
                                                      bool plasticity);
} // namespace seissol::waveprop
