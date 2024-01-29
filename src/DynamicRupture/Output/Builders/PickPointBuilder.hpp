#ifndef SEISSOL_DR_OUTPUT_PICKPOINT_BUILDER_HPP
#define SEISSOL_DR_OUTPUT_PICKPOINT_BUILDER_HPP

#include "ReceiverBasedOutputBuilder.hpp"
#include <Initializer/PointMapper.h>

namespace seissol::dr::output {
class PickPointBuilder : public ReceiverBasedOutputBuilder {
  public:
  ~PickPointBuilder() override = default;
  void setParams(PickpointParams params) { pickpointParams = std::move(params); }
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
    using namespace seissol::initializers;
    StringsType content = FileProcessor::getFileAsStrings(pickpointParams.ppFileName);
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

    // findMeshIds expects a vector of eigenPoints.
    // Therefore, we need to convert
    std::vector<Eigen::Vector3d> eigenPoints(numReceiverPoints);
    for (size_t receiverId{0}; receiverId < numReceiverPoints; ++receiverId) {
      const auto& receiverPoint = potentialReceivers[receiverId];
      eigenPoints[receiverId] = receiverPoint.global.getAsEigen3LibVector();
    }

    std::vector<short> contained(numReceiverPoints);
    std::vector<unsigned> localIds(numReceiverPoints, std::numeric_limits<unsigned>::max());

    const auto [faultVertices, faultElements, elementToFault] = getElementsAlongFault();

    seissol::initializers::findMeshIds(eigenPoints.data(),
                                       faultVertices,
                                       faultElements,
                                       numReceiverPoints,
                                       contained.data(),
                                       localIds.data());

    const auto& meshElements = meshReader->getElements();
    const auto& meshVertices = meshReader->getVertices();
    const auto& faultInfos = meshReader->getFault();

    for (size_t receiverIdx{0}; receiverIdx < numReceiverPoints; ++receiverIdx) {
      auto& receiver = potentialReceivers[receiverIdx];

      if (static_cast<bool>(contained[receiverIdx])) {
        const auto localId = localIds[receiverIdx];

        const auto faultIndicesIt = elementToFault.find(localId);
        assert(faultIndicesIt != elementToFault.end());
        const auto& faultIndices = faultIndicesIt->second;

        const auto firstFaultIdx = faultIndices[0];

        // find the original element which contains a fault face
        // note: this allows to project a receiver to the plus side
        //       even if it was found in the negative one
        const auto element = meshElements[faultInfos[firstFaultIdx].element];

        const auto closest = findClosestFaultIndex(receiver.global, element, faultIndices);
        const auto faultItem = faultInfos[closest];

        receiver.globalTriangle = getGlobalTriangle(faultItem.side, element, meshVertices);
        projectPointToFace(receiver.global, receiver.globalTriangle, faultItem.normal);

        receiver.isInside = true;
        receiver.faultFaceIndex = closest;
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

  std::tuple<std::vector<Vertex>,
             std::vector<Element>,
             std ::unordered_map<size_t, std::vector<size_t>>>
      getElementsAlongFault() {
    const auto& fault = meshReader->getFault();
    const auto numFaultElements = fault.size();

    auto meshElements = meshReader->getElements();
    auto meshVertices = meshReader->getVertices();

    constexpr int numSides{2};
    constexpr int numVertices{4};

    std::vector<Vertex> faultVertices;
    faultVertices.reserve(numVertices * numFaultElements * numSides);

    std::vector<Element> faultElements;
    faultElements.reserve(numFaultElements * numSides);

    std::vector<size_t> filtertedToOrigIndices;
    filtertedToOrigIndices.reserve(2 * numFaultElements);

    // note: an element can have multiple fault faces
    std::unordered_map<size_t, std::vector<size_t>> elementToFault{};

    for (size_t faultIdx{0}; faultIdx < numFaultElements; ++faultIdx) {
      const auto& faultItem = fault[faultIdx];

      if ((faultItem.element >= 0) && (faultItem.neighborElement >= 0)) {
        // element copy done on purpose because we are recording
        // a new vertex array and thus we need to modify vertex indices
        // inside of each element
        std::array<Element*, numSides> elements{&meshElements[faultItem.element],
                                                &meshElements[faultItem.neighborElement]};
        for (int sideId = 0; sideId < numSides; ++sideId) {
          auto element = (*elements[sideId]);

          for (size_t vertexIdx{0}; vertexIdx < numVertices; ++vertexIdx) {
            faultVertices.push_back(meshVertices[element.vertices[vertexIdx]]);
            element.vertices[vertexIdx] = faultVertices.size() - 1;
          }
          faultElements.push_back(element);
          elementToFault[element.localId].push_back(faultIdx);
        }
      }
    }
    faultVertices.shrink_to_fit();
    faultElements.shrink_to_fit();
    return std::make_tuple(faultVertices, faultElements, elementToFault);
  }

  size_t findClosestFaultIndex(const ExtVrtxCoords& point,
                               const Element& element,
                               const std::vector<size_t>& faultIndices) {
    assert(!faultIndices.empty() && "an element must contain some rupture faces");

    // Note: it is not so common to have an element with multiple rupture faces.
    //       Handling a trivial solution
    if (faultIndices.size() == 1) {
      return faultIndices[0];
    }

    const auto meshVertices = meshReader->getVertices();
    const auto& fault = meshReader->getFault();

    auto minDistance = std::numeric_limits<double>::max();
    auto closest = std::numeric_limits<size_t>::max();

    for (auto faceIdx : faultIndices) {
      const auto& faultItem = fault[faceIdx];

      const auto face = getGlobalTriangle(faultItem.side, element, meshVertices);

      const auto distance = getDistanceFromPointToFace(point, face, faultItem.normal);
      if (minDistance > distance) {
        minDistance = distance;
        closest = faceIdx;
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
      for (size_t idx{0}; idx < size; ++idx) {
        const auto isFound = globalContainVector[idx];
        if (!isFound) {
          logInfo(localRank) << "pickpoint fault output: "
                             << "receiver (" << idx + 1 << ") is not inside "
                             << "any element along the rupture surface";
          allReceiversFound = false;
        }
      }
      if (allReceiversFound) {
        logInfo(localRank) << "all point receivers found along the fault";
      }
    }
  }

  private:
  PickpointParams pickpointParams;
  std::vector<ReceiverPoint> potentialReceivers{};
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_PICKPOINT_BUILDER_HPP
