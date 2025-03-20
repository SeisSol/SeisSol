// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_PICKPOINTBUILDER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_PICKPOINTBUILDER_H_

#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/PointMapper.h"
#include "ReceiverBasedOutputBuilder.h"

#include <Common/Iterator.h>
#include <optional>

namespace seissol::dr::output {
class PickPointBuilder : public ReceiverBasedOutputBuilder {
  public:
  ~PickPointBuilder() override = default;
  void setParams(seissol::initializer::parameters::PickpointParameters params) {
    pickpointParams = std::move(params);
  }
  void build(std::shared_ptr<ReceiverOutputData> pickPointOutputData) override {
    outputData = pickPointOutputData;
    readCoordsFromFile();
    initReceiverLocations();
    assignNearestGaussianPoints(outputData->receiverPoints);
    assignNearestInternalGaussianPoints();
    assignFaultTags();
    initTimeCaching();
    initFaultDirections();
    initOutputVariables(pickpointParams.outputMask);
    initRotationMatrices();
    initBasisFunctions();
    initJacobian2dMatrices();
    outputData->isActive = true;
  }

  protected:
  void readCoordsFromFile() {
    using namespace seissol::initializer;
    StringsType content = FileProcessor::getFileAsStrings(pickpointParams.pickpointFileName,
                                                          "pickpoint/on-fault receiver file");
    FileProcessor::removeEmptyLines(content);

    // iterate line by line and initialize DrRecordPoints
    for (const auto& line : content) {
      std::array<real, 3> coords{};
      convertStringToMask(line, coords);

      ReceiverPoint point{};
      for (int i = 0; i < 3; ++i) {
        point.global.coords[i] = coords[i];
      }

      potentialReceivers.push_back(point);
    }
  }

  void initReceiverLocations() {
    const auto numReceiverPoints = potentialReceivers.size();

    const auto& meshElements = meshReader->getElements();
    const auto& meshVertices = meshReader->getVertices();
    const auto& faultInfos = meshReader->getFault();

    std::vector<short> contained(potentialReceivers.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t receiverIdx = 0; receiverIdx < numReceiverPoints; ++receiverIdx) {
      auto& receiver = potentialReceivers[receiverIdx];

      const auto closest = findClosestFaultIndex(receiver.global);

      if (closest.has_value()) {
        const auto& faultItem = faultInfos.at(closest.value());
        const auto& element = meshElements.at(faultItem.element);

        receiver.globalTriangle = getGlobalTriangle(faultItem.side, element, meshVertices);
        projectPointToFace(receiver.global, receiver.globalTriangle, faultItem.normal);

        contained[receiverIdx] = true;
        receiver.isInside = true;
        receiver.faultFaceIndex = closest.value();
        receiver.localFaceSideId = faultItem.side;
        receiver.globalReceiverIndex = receiverIdx;
        receiver.elementIndex = element.localId;

        receiver.reference =
            transformations::tetrahedronGlobalToReference(meshVertices[element.vertices[0]].coords,
                                                          meshVertices[element.vertices[1]].coords,
                                                          meshVertices[element.vertices[2]].coords,
                                                          meshVertices[element.vertices[3]].coords,
                                                          receiver.global.getAsEigen3LibVector());
      }
    }

    reportFoundReceivers(contained);
    for (auto& receiver : potentialReceivers) {
      if (receiver.isInside) {
        outputData->receiverPoints.push_back(receiver);
      }
    }
  }

  std::optional<size_t> findClosestFaultIndex(const ExtVrtxCoords& point) {
    const auto& meshElements = meshReader->getElements();
    const auto& meshVertices = meshReader->getVertices();
    const auto& fault = meshReader->getFault();

    auto minDistance = std::numeric_limits<double>::max();
    auto closest = std::optional<std::size_t>();

    for (auto [faceIdx, faultItem] : seissol::common::enumerate(fault)) {
      const auto face =
          getGlobalTriangle(faultItem.side, meshElements.at(faultItem.element), meshVertices);
      const auto insideQuantifier = isInsideFace(point, face, faultItem.normal);

      if (insideQuantifier > -1e-12) {
        const auto distance = getDistanceFromPointToFace(point, face, faultItem.normal);
        if (minDistance > distance) {
          minDistance = distance;
          closest = faceIdx;
        }
      }
    }
    return closest;
  }

  void initTimeCaching() override {
    outputData->maxCacheLevel = pickpointParams.maxPickStore;
    outputData->cachedTime.resize(outputData->maxCacheLevel, 0.0);
    outputData->currentCacheLevel = 0;
  }

  protected:
  void reportFoundReceivers(std::vector<short>& localContainVector) {
    const auto size = localContainVector.size();
    std::vector<short> globalContainVector(size);

    auto comm = MPI::mpi.comm();
    MPI_Reduce(const_cast<short*>(&localContainVector[0]),
               const_cast<short*>(&globalContainVector[0]),
               size,
               MPI_SHORT,
               MPI_SUM,
               0,
               comm);

    if (localRank == 0) {
      bool allReceiversFound{true};
      std::size_t missing = 0;
      for (size_t idx{0}; idx < size; ++idx) {
        const auto isFound = globalContainVector[idx];
        if (!isFound) {
          logWarning() << "On-fault receiver " << idx
                       << " is not inside any element along the rupture surface";
          allReceiversFound = false;
          ++missing;
        }
      }
      if (allReceiversFound) {
        logInfo() << "All point receivers found along the fault";
      } else {
        logError() << missing << "on-fault receivers have not been found.";
      }
    }
  }

  private:
  seissol::initializer::parameters::PickpointParameters pickpointParams;
  std::vector<ReceiverPoint> potentialReceivers{};
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_PICKPOINTBUILDER_H_
