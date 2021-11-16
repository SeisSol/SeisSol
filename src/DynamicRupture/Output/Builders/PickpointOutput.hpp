#ifndef SEISSOL_PICKPOINTOUTPUT_HPP
#define SEISSOL_PICKPOINTOUTPUT_HPP

#include "DynamicRupture/Output/DataTypes.hpp"
#include "Geometry/MeshReader.h"
#include <Initializer/InputAux.hpp>


namespace seissol::dr::output {
class PickpointOutput {
  public:
  friend Base;

  ~PickpointOutput() {
    auto deallocateVars = [](auto& var, int) { var.releaseData(); };
    aux::forEach(outputData.vars, deallocateVars);
  }

  void setParams(const PickpointParamsT& params, const MeshReader* reader) {
    pickpointParams = params;
    meshReader = reader;
  }

  void init() {
    readCoordsFromFile();
    initReceiverLocations();
    initFaultDirections();
    initTimeCaching();
    initOutputVariables();
    // findElementsContainingPoints();
    // initPointsIndices();
    // projectPointsToFaces();
    // findClosestGpPoint();
    outputData.isActive = true;
  }

  void readCoordsFromFile() {
    using namespace seissol::initializers;
    StringsT content = FileProcessor::getFileAsStrings(pickpointParams.ppFileName);
    FileProcessor::removeEmptyLines(content);

    if (pickpointParams.numOutputPoints > static_cast<int>(content.size()))
      throw std::runtime_error(
          "requested num. of fault pick-points is more than the file contains");

    // iterate line by line and initialize DrRecordPoints
    for (const auto& line : content) {
      std::array<real, 3> coords{};
      convertStringToMask(line, coords);

      ReceiverPointT point{};
      for (int i = 0; i < 3; ++i) {
        point.global.coords[i] = coords[i];
      }

      outputData.receiverPoints.push_back(point);
    }
  }

  void initReceiverLocations() {
    for (size_t i = 0; i < outputData.receiverPoints.size(); ++i) {
      outputData.receiverPoints[i].globalReceiverIndex = i;
      outputData.receiverPoints[i].faultFaceIndex = 0;
    }
  }

  void initTimeCaching() {
    outputData.maxCacheLevel = pickpointParams.maxPickStore;
    outputData.cachedTime.resize(outputData.maxCacheLevel, 0.0);
    outputData.currentCacheLevel = 0;
  }

  void initOutputVariables() {

    auto assignMask = [this](auto& var, int index) {
      var.isActive = this->pickpointParams.outputMask[index];
    };
    aux::forEach(outputData.vars, assignMask);

    auto allocateVariables = [this](auto& var, int) {
      var.maxCacheLevel = outputData.maxCacheLevel;
      var.allocateData(this->outputData.receiverPoints.size());
    };
    aux::forEach(outputData.vars, allocateVariables);
  }

  void initFaultDirections() {
    size_t size = outputData.receiverPoints.size();
    outputData.faultDirections.resize(outputData.receiverPoints.size());
    const auto& faultInfo = meshReader->getFault();

    for (size_t index = 0; index < size; ++index) {
      size_t globalIndex = outputData.receiverPoints[index].faultFaceIndex;

      outputData.faultDirections[index].faceNormal = faultInfo[globalIndex].normal;
      outputData.faultDirections[index].tangent1 = faultInfo[globalIndex].tangent1;
      outputData.faultDirections[index].tangent2 = faultInfo[globalIndex].tangent2;
      computeStrikeAndDipVectors(outputData.faultDirections[index].faceNormal,
                                 outputData.faultDirections[index].strike,
                                 outputData.faultDirections[index].dip);
    }
  }

  private:
  PickpointParamsT pickpointParams;
  const MeshReader* meshReader;
  OutputData outputData;
};
}
#endif //SEISSOL_PICKPOINTOUTPUT_HPP
