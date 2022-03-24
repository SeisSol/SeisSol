#include "DynamicRupture/Output/IntegratedOutput.hpp"

namespace seissol::dr::output {
void IntegratedOutput::setLtsData(seissol::initializers::DynamicRupture* userDrDescr,
                                  size_t numElements) {
  drDescr = userDrDescr;
  numFaultElements = numElements;
}

double IntegratedOutput::getMagnitude(GeoOutputData& outputData) {
  double magnitude{};
  for (size_t faceIndex = 0; faceIndex < numFaultElements; ++faceIndex) {
    auto ltsMap = (*faceToLtsMap)[faceIndex];
    auto* layer = ltsMap.first;
    auto ltsId = ltsMap.second;

    auto averagedSlip = (layer->var(drDescr->averagedSlip))[ltsId];
    auto* mu = (layer->var(drDescr->mu))[ltsId];

    double averageMu{0.0};
    for (size_t point = 0; point < misc::numberOfBoundaryGaussPoints; ++point) {
      averageMu += mu[point];
    }

    magnitude += averagedSlip * averageMu * outputData.surfaceAreas[faceIndex];
  }

  return magnitude;
}

double IntegratedOutput::getMomentRate(GeoOutputData& outputData) {

  real momentRate{0.0};
  for (size_t faceIndex = 0; faceIndex < numFaultElements; ++faceIndex) {
    auto ltsMap = (*faceToLtsMap)[faceIndex];
    auto* layer = ltsMap.first;
    auto ltsId = ltsMap.second;

    auto* slipRate1 = (layer->var(drDescr->slipRate1))[ltsId];
    auto* slipRate2 = (layer->var(drDescr->slipRate2))[ltsId];
    auto* mu = (layer->var(drDescr->mu))[ltsId];

    double averageSr{0.0};
    double averageMu{0.0};

    for (size_t point = 0; point < misc::numberOfBoundaryGaussPoints; ++point) {
      averageMu += mu[point];

      real slip = slipRate1[point] * slipRate1[point] + slipRate2[point] * slipRate2[point];
      averageSr += std::sqrt(slip) / static_cast<double>(misc::numberOfBoundaryGaussPoints);
    }
    averageMu /= static_cast<double>(misc::numberOfBoundaryGaussPoints);
    momentRate += averageSr * averageMu * outputData.surfaceAreas[faceIndex];
  }
  return momentRate;
}
} // namespace seissol::dr::output