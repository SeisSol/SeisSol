#ifndef SEISSOL_DR_OUTPUT_ELEMENTWISE_BUILDER_HPP
#define SEISSOL_DR_OUTPUT_ELEMENTWISE_BUILDER_HPP

#include "DynamicRupture/Output/FaultRefiner/FaultRefiners.hpp"
#include "OutputBuilder.hpp"

namespace seissol::dr::output {
class Base;

class ElementWiseBuilder : public OutputBuilder {
  public:
  friend Base;

  ~ElementWiseBuilder() override = default;

  void setParams(const ElementwiseFaultParamsT& params) { elementwiseParams = params; }

  void init() override {
    initReceiverLocations();
    assignNearestGaussianPoints(outputData.receiverPoints);
    initTimeCaching();
    initOutputVariables(this->elementwiseParams.outputMask);
    initFaultDirections();
    initRotationMatrices();
    initBasisFunctions();
    outputData.isActive = true;
  }

  void initTimeCaching() override {
    outputData.maxCacheLevel = ElementWiseBuilder::maxAllowedCacheLevel;
    outputData.currentCacheLevel = 0;
  }

  void initReceiverLocations() {
    std::unique_ptr<refiner::FaultRefiner> faultRefiner{nullptr};
    faultRefiner = refiner::get(elementwiseParams.refinementStrategy);

    const auto numFaultElements = meshReader->getFault().size();
    const auto numSubTriangles = faultRefiner->getNumSubTriangles();

    logInfo(localRank) << "CPP: Initialising Fault output. "
                       << "Number of sub-triangles: " << numSubTriangles;

    // get arrays of elements and vertices from the mesher
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    // iterate through each fault side
    for (size_t faceIndex = 0; faceIndex < numFaultElements; ++faceIndex) {

      // get a Global Element ID for the current fault face
      auto elementIndex = faultInfo[faceIndex].element;
      const auto& element = elementsInfo[elementIndex];

      // store coords of vertices of the current ELEMENT
      std::array<const double*, 4> elementVerticesCoords{};
      for (int elementVertexId = 0; elementVertexId < 4; ++elementVertexId) {
        auto globalVertexId = element.vertices[elementVertexId];
        elementVerticesCoords[elementVertexId] = verticesInfo[globalVertexId].coords;
      }

      auto localFaceSideId = faultInfo[faceIndex].side;

      // init reference coordinates of the fault face
      ExtTriangle referenceFace = getReferenceFace(localFaceSideId);

      // init global coordinates of the fault face
      ExtTriangle globalFace = getGlobalTriangle(localFaceSideId, element, verticesInfo);

      faultRefiner->refineAndAccumulate(
          {elementwiseParams.refinement, static_cast<int>(faceIndex), localFaceSideId},
          std::make_pair(globalFace, referenceFace));
    }

    // retrieve all receivers from a fault face refiner
    outputData.receiverPoints = faultRefiner->moveAllReceiverPoints();
    faultRefiner.reset(nullptr);
  }

  void initConstrains() {}

  inline const static size_t maxAllowedCacheLevel = 1;

  private:
  ElementwiseFaultParamsT elementwiseParams;
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_ELEMENTWISE_BUILDER_HPP
