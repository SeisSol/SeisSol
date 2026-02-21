// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_ELEMENTWISEBUILDER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_ELEMENTWISEBUILDER_H_

#include "DynamicRupture/Output/FaultRefiner/FaultRefiners.h"
#include "DynamicRupture/Output/Geometry.h"
#include "DynamicRupture/Output/OutputAux.h"
#include "GeneratedCode/init.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Numerical/Transformation.h"
#include "ReceiverBasedOutputBuilder.h"
#include "Solver/MultipleSimulations.h"

namespace seissol::dr::output {
class ElementWiseBuilder : public ReceiverBasedOutputBuilder {
  public:
  ~ElementWiseBuilder() override = default;
  void setParams(const seissol::initializer::parameters::ElementwiseFaultParameters& params) {
    elementwiseParams = params;
  }
  void build(std::shared_ptr<ReceiverOutputData> elementwiseOutputData) {
    outputData = std::move(elementwiseOutputData);
    initReceiverLocations();
    assignNearestGaussianPoints(outputData->receiverPoints);
    assignNearestInternalGaussianPoints();
    assignFusedIndices();
    assignFaultTags();
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
    outputData->maxCacheLevel = ElementWiseBuilder::MaxAllowedCacheLevel;
    outputData->currentCacheLevel = 0;
  }

  void initReceiverLocations() {
    auto faultRefiner = refiner::get(elementwiseParams.refinementStrategy);

    const auto numSubTriangles = faultRefiner->getNumSubTriangles();
    const auto order = static_cast<std::uint32_t>(std::max(elementwiseParams.vtkorder, 0));

    const auto numFaultElements = meshReader->getFault().size();

    logInfo() << "Initializing Fault output."
              << "Number of sub-triangles:" << numSubTriangles << "Output order:" << order
              << "Simulation count:" << multisim::NumSimulations;

    // get arrays of elements and vertices from the meshReader
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    // iterate through each fault side
    for (size_t faceIdx = 0; faceIdx < numFaultElements; ++faceIdx) {
      const auto& fault = faultInfo[faceIdx];
      const auto elementIdx = fault.element;

      if (elementIdx >= 0) {
        const auto& element = elementsInfo[elementIdx];

        // store coords of vertices of the current ELEMENT
        std::array<const double*, Cell::NumVertices> elementVerticesCoords{};
        for (size_t vertexIdx = 0; vertexIdx < Cell::NumVertices; ++vertexIdx) {
          auto globalVertexIdx = element.vertices[vertexIdx];
          elementVerticesCoords[vertexIdx] = verticesInfo[globalVertexIdx].coords;
        }

        const auto faceSideIdx = fault.side;

        // init reference coordinates of the fault face
        const ExtTriangle referenceTriangle = getReferenceTriangle(faceSideIdx);

        // init global coordinates of the fault face
        const ExtTriangle globalFace = getGlobalTriangle(faceSideIdx, element, verticesInfo);

        faultRefiner->refineAndAccumulate({elementwiseParams.refinement,
                                           static_cast<int>(faceIdx),
                                           faceSideIdx,
                                           elementIdx,
                                           element.globalId,
                                           order,
                                           multisim::NumSimulations},
                                          std::make_pair(globalFace, referenceTriangle));
      }
    }

    // retrieve all receivers from a fault face refiner
    outputData->receiverPoints = faultRefiner->moveAllReceiverPoints();
    faultRefiner.reset(nullptr);
  }

  inline const static size_t MaxAllowedCacheLevel = 1;

  private:
  seissol::initializer::parameters::ElementwiseFaultParameters elementwiseParams;
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_ELEMENTWISEBUILDER_H_
