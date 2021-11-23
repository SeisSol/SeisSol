#ifndef SEISSOL_DR_PICKPOINT_BUILDER_HPP
#define SEISSOL_DR_PICKPOINT_BUILDER_HPP

#include "OutputBuilder.hpp"
#include <Initializer/PointMapper.h>

namespace seissol::dr::output {
class PickPointBuilder : public OutputBuilder {
  public:
  friend Base;

  ~PickPointBuilder() override = default;

  void setParams(const PickpointParamsT& params) { pickpointParams = params; }

  void init() override {
    readCoordsFromFile();
    initReceiverLocations();
    assignNearestGaussianPoints(outputData.receiverPoints);
    initTimeCaching();
    initFaultDirections();
    initOutputVariables(pickpointParams.outputMask);
    initRotationMatrices();
    initBasisFunctions();

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
    // TODO: split method into smaller ones

    const auto numPoints = outputData.receiverPoints.size();

    std::vector<Eigen::Vector3d> eigenPoints(numPoints);
    for (size_t pointIdx{0}; pointIdx < numPoints; ++pointIdx) {
      auto& receiverPoint = outputData.receiverPoints[pointIdx];
      eigenPoints[pointIdx] = receiverPoint.global.getAsEigenVector();
    }

    std::vector<short> contained(numPoints);
    std::vector<unsigned> meshIds(numPoints, std::numeric_limits<unsigned>::max());

    std::vector<Vertex> faultVertices;
    std::vector<Element> faultElements;
    std::unordered_map<size_t, std::vector<size_t>> elementToFault{};

    std::tie(faultVertices, faultElements, elementToFault) = getElementsAlongFault();

    seissol::initializers::findMeshIds(eigenPoints.data(),
                                       faultVertices,
                                       faultElements,
                                       numPoints,
                                       contained.data(),
                                       meshIds.data());

    auto meshElements = meshReader->getElements();
    auto meshVertices = meshReader->getVertices();

    const auto& fault = meshReader->getFault();
    for (size_t receiverIdx{0}; receiverIdx < numPoints; ++receiverIdx) {
      auto& receiver = outputData.receiverPoints[receiverIdx];

      if (static_cast<bool>(contained[0])) {
        auto meshId = meshIds[receiverIdx];
        const auto& faultIndices = elementToFault[meshId];

        // TODO: needs to be extended to a general case where an element can have
        //       more than one face rupture
        assert(faultIndices.size() == 1 &&
               "currently, we assume that only one element contains a DR receiver");
        auto faultItem = fault[faultIndices[0]];

        receiver.globalSubTet = getGlobalFace(faultItem, meshElements, meshVertices);
        projectPointToFace(receiver.global, receiver.globalSubTet, faultItem.normal);

        receiver.isInside = true;
        receiver.faultFaceIndex = faultIndices[0];
        receiver.localFaceSideId = faultItem.side;
        receiver.globalReceiverIndex = receiverIdx;
        receiver.elementIndex = meshId;

        using namespace seissol::transformations;
        auto element = meshElements[meshId];
        receiver.referece = tetrahedronGlobalToReference(meshVertices[element.vertices[0]].coords,
                                                         meshVertices[element.vertices[1]].coords,
                                                         meshVertices[element.vertices[2]].coords,
                                                         meshVertices[element.vertices[3]].coords,
                                                         receiver.global.getAsEigenVector());

        receiver.isInside = true;

      } else {
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

    std::vector<Vertex> faultVertices(4 * numFaultElements);
    std::vector<Element> faultElements(numFaultElements);

    // note: an element can have multiple fault faces
    std::unordered_map<size_t, std::vector<size_t>> elementToFault{};

    for (size_t faultIdx{0}; faultIdx < numFaultElements; ++faultIdx) {
      const auto& faultItem = fault[faultIdx];
      Element element = meshElements[faultItem.element];

      for (size_t vertexIdx{0}; vertexIdx < 4; ++vertexIdx) {
        Vertex vertex = meshVertices[element.vertices[vertexIdx]];
        const size_t faultVertexIdx = vertexIdx + 4 * faultIdx;
        faultVertices[faultVertexIdx] = vertex;

        // note: we need to remap vertex indicies inside of an element
        // because we record a new vertex vector
        element.vertices[vertexIdx] = faultVertexIdx;
      }

      faultElements[faultIdx] = element;
      elementToFault[faultItem.element].push_back(faultIdx);
    }

    return std::make_tuple(faultVertices, faultElements, elementToFault);
  }

  void initTimeCaching() override {
    outputData.maxCacheLevel = pickpointParams.maxPickStore;
    outputData.cachedTime.resize(outputData.maxCacheLevel, 0.0);
    outputData.currentCacheLevel = 0;
  }

  private:
  PickpointParamsT pickpointParams;
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_PICKPOINT_BUILDER_HPP
