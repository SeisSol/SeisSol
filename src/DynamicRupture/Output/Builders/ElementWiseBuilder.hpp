#ifndef SEISSOL_DR_OUTPUT_ELEMENTWISE_BUILDER_HPP
#define SEISSOL_DR_OUTPUT_ELEMENTWISE_BUILDER_HPP

#include "DynamicRupture/Output/FaultRefiner/FaultRefiners.hpp"
#include "ReceiverBasedOutputBuilder.hpp"

namespace seissol::dr::output {
class ElementWiseBuilder : public ReceiverBasedOutputBuilder {
  public:
  ~ElementWiseBuilder() override = default;
  void setParams(const ElementwiseFaultParamsT& params) { elementwiseParams = params; }
  void build(ReceiverBasedOutputData* ewOutputData) override {
    outputData = ewOutputData;
    initReceiverLocations();
    assignNearestGaussianPoints(outputData->receiverPoints);
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
    std::unique_ptr<refiner::FaultRefiner> faultRefiner{nullptr};
    faultRefiner = refiner::get(elementwiseParams.refinementStrategy);

    const auto numFaultElements = meshReader->getFault().size();
    const auto numSubTriangles = faultRefiner->getNumSubTriangles();

    logInfo(localRank) << "CPP: Initialising Fault output. "
                       << "Number of sub-triangles: " << numSubTriangles;

    // get arrays of elements and vertices from the meshReader
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    // iterate through each fault side
    for (size_t faceIndex = 0; faceIndex < numFaultElements; ++faceIndex) {

      // get a global element ID for the current fault face
      const auto& fault = faultInfo[faceIndex];
      auto elementIndex = fault.element;

      if (elementIndex >= 0) {
        const auto& element = elementsInfo[elementIndex];

        // store coords of vertices of the current ELEMENT
        constexpr size_t numSides{4};
        std::array<const double*, numSides> elementVerticesCoords{};
        for (size_t side = 0; side < numSides; ++side) {
          auto globalVertexId = element.vertices[side];
          elementVerticesCoords[side] = verticesInfo[globalVertexId].coords;
        }

        auto faceSideId = fault.side;

        // init reference coordinates of the fault face
        ExtTriangle referenceFace = getReferenceFace(faceSideId);

        // init global coordinates of the fault face
        ExtTriangle globalFace = getGlobalTriangle(faceSideId, element, verticesInfo);

        faultRefiner->refineAndAccumulate(
            {elementwiseParams.refinement, static_cast<int>(faceIndex), faceSideId, elementIndex},
            std::make_pair(globalFace, referenceFace));
      }
    }

    // retrieve all receivers from a fault face refiner
    outputData->receiverPoints = faultRefiner->moveAllReceiverPoints();
    faultRefiner.reset(nullptr);
  }

  inline const static size_t maxAllowedCacheLevel = 1;

  private:
  ElementwiseFaultParamsT elementwiseParams;
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_ELEMENTWISE_BUILDER_HPP
