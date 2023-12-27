#ifndef SEISSOL_DR_OUTPUT_ELEMENTWISE_BUILDER_HPP
#define SEISSOL_DR_OUTPUT_ELEMENTWISE_BUILDER_HPP

#include "DynamicRupture/Output/FaultRefiner/FaultRefiners.hpp"
#include "ReceiverBasedOutputBuilder.hpp"
#include "Initializer/parameters/OutputParameters.h"

namespace seissol::dr::output {
class ElementWiseBuilder : public ReceiverBasedOutputBuilder {
  public:
  ~ElementWiseBuilder() override = default;
  void setParams(const seissol::initializers::parameters::ElementwiseFaultParameters& params) {
    elementwiseParams = params;
  }
  void build(std::shared_ptr<ReceiverOutputData> elementwiseOutputData) override {
    outputData = elementwiseOutputData;
    initReceiverLocations();
    assignNearestGaussianPoints(outputData->receiverPoints);
    assignNearestInternalGaussianPoints();
    initTimeCaching();
    initOutputVariables(elementwiseParams.outputMask);
    initFaultDirections();
    initRotationMatrices();
    initBasisFunctions();
    initJacobian2dMatrices();
    outputData->isActive = true;
  }

  protected:
  void initTimeCaching() override {
    outputData->maxCacheLevel = ElementWiseBuilder::maxAllowedCacheLevel;
    outputData->currentCacheLevel = 0;
  }

  void initReceiverLocations() {
    auto faultRefiner = refiner::get(elementwiseParams.refinementStrategy);

    const auto numFaultElements = meshReader->getFault().size();
    const auto numSubTriangles = faultRefiner->getNumSubTriangles();

    logInfo(localRank) << "Initializing Fault output."
                       << "Number of sub-triangles:" << numSubTriangles;

    // get arrays of elements and vertices from the meshReader
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    // iterate through each fault side
    for (size_t faceIdx = 0; faceIdx < numFaultElements; ++faceIdx) {

      // get a global element ID for the current fault face
      const auto& fault = faultInfo[faceIdx];
      const auto elementIdx = fault.element;

      if (elementIdx >= 0) {
        const auto& element = elementsInfo[elementIdx];

        // store coords of vertices of the current ELEMENT
        constexpr size_t numVertices{4};
        std::array<const double*, numVertices> elementVerticesCoords{};
        for (size_t vertexIdx = 0; vertexIdx < numVertices; ++vertexIdx) {
          auto globalVertexIdx = element.vertices[vertexIdx];
          elementVerticesCoords[vertexIdx] = verticesInfo[globalVertexIdx].coords;
        }

        const auto faceSideIdx = fault.side;

        // init reference coordinates of the fault face
        ExtTriangle referenceTriangle = getReferenceTriangle(faceSideIdx);

        // init global coordinates of the fault face
        ExtTriangle globalFace = getGlobalTriangle(faceSideIdx, element, verticesInfo);

        faultRefiner->refineAndAccumulate(
            {elementwiseParams.refinement, static_cast<int>(faceIdx), faceSideIdx, elementIdx},
            std::make_pair(globalFace, referenceTriangle));
      }
    }

    // retrieve all receivers from a fault face refiner
    outputData->receiverPoints = faultRefiner->moveAllReceiverPoints();
    faultRefiner.reset(nullptr);
  }

  inline const static size_t maxAllowedCacheLevel = 1;

  private:
  seissol::initializers::parameters::ElementwiseFaultParameters elementwiseParams;
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_ELEMENTWISE_BUILDER_HPP
