#include "DynamicRupture/Output/IntegratedOutput.hpp"

namespace seissol::dr::output {
void IntegratedOutput::setLtsData(seissol::initializers::DynamicRupture* userDrDescr,
                                  size_t numElements) {
  drDescr = userDrDescr;
  numFaultElements = numElements;
}

long double IntegratedOutput::getSeismicMoment(IntegratedOutputData& outputData) {
  long double magnitude{};
  for (size_t faceIndex = 0; faceIndex < numFaultElements; ++faceIndex) {
    auto ltsMap = (*faceToLtsMap)[faceIndex];
    auto* layer = ltsMap.first;
    auto ltsId = ltsMap.second;

    auto averagedSlip = (layer->var(drDescr->averagedSlip))[ltsId];
    auto lambda = outputData.lambda[faceIndex];
    magnitude += averagedSlip * lambda * outputData.surfaceAreas[faceIndex];
  }
  return magnitude;
}

double IntegratedOutput::getSeismicMomentRate(IntegratedOutputData& outputData) {
  double momentRate{0.0};
  for (size_t faceIndex = 0; faceIndex < numFaultElements; ++faceIndex) {
    auto ltsMap = (*faceToLtsMap)[faceIndex];
    auto* layer = ltsMap.first;
    auto ltsId = ltsMap.second;

    auto* slipRate1 = (layer->var(drDescr->slipRate1))[ltsId];
    auto* slipRate2 = (layer->var(drDescr->slipRate2))[ltsId];

    double averageSlipRate{0.0};
    for (size_t point = 0; point < misc::numberOfBoundaryGaussPoints; ++point) {
      real normSlipRateSquared =
          slipRate1[point] * slipRate1[point] + slipRate2[point] * slipRate2[point];
      averageSlipRate +=
          std::sqrt(normSlipRateSquared) / static_cast<double>(misc::numberOfBoundaryGaussPoints);
    }

    auto lambda = outputData.lambda[faceIndex];
    momentRate += averageSlipRate * lambda * outputData.surfaceAreas[faceIndex];
  }
  return momentRate;
}
} // namespace seissol::dr::output