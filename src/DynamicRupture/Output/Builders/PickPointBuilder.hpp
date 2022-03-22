#ifndef SEISSOL_DR_OUTPUT_PICKPOINT_BUILDER_HPP
#define SEISSOL_DR_OUTPUT_PICKPOINT_BUILDER_HPP

#include "OutputBuilder.hpp"
#include <Initializer/PointMapper.h>

namespace seissol::dr::output {
class PickPointBuilder : public OutputBuilder {
  public:
  ~PickPointBuilder() override = default;
  void setParams(const PickpointParamsT& params) { pickpointParams = params; }
  void build(OutputData* ppOutputData) override {
    outputData = ppOutputData;
    readCoordsFromFile();
    initReceiverLocations();
    assignNearestGaussianPoints(outputData->receiverPoints);
    initTimeCaching();
    initFaultDirections();
    initOutputVariables(pickpointParams.outputMask);
    initRotationMatrices();
    initBasisFunctions();
    initJacobian2dMatrices();
    outputData->isActive = true;
  }

  void readCoordsFromFile() {
    using namespace seissol::initializers;
    StringsT content = FileProcessor::getFileAsStrings(pickpointParams.ppFileName);
    FileProcessor::removeEmptyLines(content);

    // iterate line by line and initialize DrRecordPoints
    for (const auto& line : content) {
      std::array<real, 3> coords{};
      convertStringToMask(line, coords);

      ReceiverPointT point{};
      for (int i = 0; i < 3; ++i) {
        point.global.coords[i] = coords[i];
      }

      outputData->receiverPoints.push_back(point);
    }
  }

  void initReceiverLocations() {
    const auto numReceiverPoints = outputData->receiverPoints.size();

    // findMeshIds expects a vector of eigenPoints.
    // Therefore, we need to convert
    std::vector<Eigen::Vector3d> eigenPoints(numReceiverPoints);
    for (size_t receiverId{0}; receiverId < numReceiverPoints; ++receiverId) {
      auto& receiverPoint = outputData->receiverPoints[receiverId];
      eigenPoints[receiverId] = receiverPoint.global.getAsEigenVector();
    }

    std::vector<short> contained(numReceiverPoints);
    std::vector<unsigned> meshIds(numReceiverPoints, std::numeric_limits<unsigned>::max());

    std::vector<Vertex> faultVertices;
    std::vector<Element> faultElements;
    std::unordered_map<size_t, std::vector<size_t>> elementToFault{};

    std::tie(faultVertices, faultElements, elementToFault) = getElementsAlongFault();

    seissol::initializers::findMeshIds(eigenPoints.data(),
                                       faultVertices,
                                       faultElements,
                                       numReceiverPoints,
                                       contained.data(),
                                       meshIds.data());

    auto meshElements = meshReader->getElements();
    auto meshVertices = meshReader->getVertices();

    const auto& fault = meshReader->getFault();
    for (size_t receiverId{0}; receiverId < numReceiverPoints; ++receiverId) {
      auto& receiver = outputData->receiverPoints[receiverId];

      if (static_cast<bool>(contained[receiverId])) {
        auto meshId = meshIds[receiverId];
        const auto& faultIndices = elementToFault[meshId];

        auto element = meshElements[meshId];
        auto closest = findClosestFaultIndex(receiver.global, element, faultIndices);
        auto faultItem = fault[closest];

        receiver.globalTriangle = getGlobalTriangle(faultItem.side, element, meshVertices);
        projectPointToFace(receiver.global, receiver.globalTriangle, faultItem.normal);

        receiver.isInside = true;
        receiver.faultFaceIndex = closest;
        receiver.localFaceSideId = faultItem.side;
        receiver.globalReceiverIndex = receiverId;
        receiver.elementIndex = meshId;

        using namespace seissol::transformations;
        receiver.reference = tetrahedronGlobalToReference(meshVertices[element.vertices[0]].coords,
                                                          meshVertices[element.vertices[1]].coords,
                                                          meshVertices[element.vertices[2]].coords,
                                                          meshVertices[element.vertices[3]].coords,
                                                          receiver.global.getAsEigenVector());
        receiver.isInside = true;
      } else {
        logInfo() << "pickpoint fault output: receiver (" << receiverId << ") is not inside"
                  << "of any element along the rupture surface";
        receiver.isInside = false;
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
    std::vector<Vertex> faultVertices(4 * numFaultElements * numSides);
    std::vector<Element> faultElements(numFaultElements * numSides);

    // note: an element can have multiple fault faces
    std::unordered_map<size_t, std::vector<size_t>> elementToFault{};

    for (size_t faultIdx{0}; faultIdx < numFaultElements; ++faultIdx) {
      const auto& faultItem = fault[faultIdx];

      // element copy done on purpose because we are recording
      // a new vertex array and thus we need to modify vertex indices
      // inside of each element
      std::array<Element*, numSides> elements{&meshElements[faultItem.element],
                                              &meshElements[faultItem.neighborElement]};
      for (int sideId = 0; sideId < numSides; ++sideId) {
        auto& element = elements[sideId];

        for (size_t vertexIdx{0}; vertexIdx < 4; ++vertexIdx) {
          const size_t faultVertexIdx = vertexIdx + 4 * (sideId + numSides * faultIdx);
          faultVertices[faultVertexIdx] = meshVertices[element->vertices[vertexIdx]];
          element->vertices[vertexIdx] = faultVertexIdx;
        }

        faultElements[sideId + numSides * faultIdx] = *element;
      }
      elementToFault[faultItem.element].push_back(faultIdx);
      elementToFault[faultItem.neighborElement].push_back(faultIdx);
    }

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

    const auto& fault = meshReader->getFault();
    const auto& meshVertices = meshReader->getVertices();

    auto minDistance = std::numeric_limits<double>::max();
    auto closest = std::numeric_limits<size_t>::max();

    for (auto faceIdx : faultIndices) {
      const auto& faultItem = fault[faceIdx];

      auto face = getGlobalTriangle(faultItem.side, element, meshVertices);

      auto distance = getDistanceFromPointToFace(point, face, faultItem.normal);
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

  private:
  PickpointParamsT pickpointParams;
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_PICKPOINT_BUILDER_HPP
