#ifndef SEISSOL_DR_ELEMENTWISE_BUILDER_HPP
#define SEISSOL_DR_ELEMENTWISE_BUILDER_HPP

#include "OutputBuilder.hpp"
#include <DynamicRupture/Math.h>
#include "DynamicRupture/Output/FaultRefiner/FaultRefiner.hpp"

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
    std::unique_ptr<FaultRefinerInterface> faultRefiner{nullptr};
    faultRefiner = getRefiner(elementwiseParams.refinementStrategy);

    const auto size = meshReader->getFault().size();
    const auto numSubTriangles = faultRefiner->getNumSubTriangles();

    logInfo(localRank) << "CPP: Initialising Fault output. Refinement strategy: "
                       << elementwiseParams.refinementStrategy
                       << " Number of sub-triangles: " << numSubTriangles;

    // get arrays of elements and vertices from the mesher
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    // iterate through each fault side
    for (size_t faceIndex = 0; faceIndex < size; ++faceIndex) {

      // get a Global Element ID for the current fault face
      auto elementIndex = faultInfo[faceIndex].element;
      const auto& element = elementsInfo[elementIndex];

      // TODO: check whether we need this if-statement
      if (elementIndex > 0) {

        // store coords of vertices of the current ELEMENT
        std::array<const double*, 4> elementVerticesCoords{};
        for (int ElementVertexId = 0; ElementVertexId < 4; ++ElementVertexId) {
          auto globalVertexId = element.vertices[ElementVertexId];
          elementVerticesCoords[ElementVertexId] = verticesInfo[globalVertexId].coords;
        }

        auto localFaceSideId = faultInfo[faceIndex].side;

        // init reference coordinates of the fault face
        ExtTriangle referenceFace = getReferenceFace(localFaceSideId);

        // init global coordinates of the fault face
        ExtTriangle globalFace = getGlobalFace(localFaceSideId, element, verticesInfo);

        faultRefiner->refineAndAccumulate(
            elementwiseParams.refinement, faceIndex, localFaceSideId, referenceFace, globalFace);
      }
    }

    // retrieve all receivers from a fault face refiner
    outputData.receiverPoints = faultRefiner->moveAllReceiverPoints();
    faultRefiner.reset(nullptr);
  }

  void initConstrains() {
    /*
    for (const auto& Point: m_ReceiverPoints) {
      const auto& RotationMatrix = m_RotationMatrices[Point.FaultFaceIndex];

    }
    */
  }

  void evaluateInitialStressInFaultCS() {
    // Compute initialStressInFaultCS
  }

  inline const static size_t maxAllowedCacheLevel = 1;

  private:
  ElementwiseFaultParamsT elementwiseParams;
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_ELEMENTWISE_BUILDER_HPP
