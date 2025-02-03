// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
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
#include "Initializer/Parameters/OutputParameters.h"
#include "Numerical/Transformation.h"
#include "ReceiverBasedOutputBuilder.h"
#include <init.h>

namespace seissol::dr::output {
class ElementWiseBuilder : public ReceiverBasedOutputBuilder {
  public:
  ~ElementWiseBuilder() override = default;
  void setParams(const seissol::initializer::parameters::ElementwiseFaultParameters& params) {
    elementwiseParams = params;
  }
  void build(std::shared_ptr<ReceiverOutputData> elementwiseOutputData) override {
    outputData = elementwiseOutputData;
    initReceiverLocations();
    assignNearestGaussianPoints(outputData->receiverPoints);
    assignNearestInternalGaussianPoints();
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
    if (elementwiseParams.vtkorder < 0) {
      auto faultRefiner = refiner::get(elementwiseParams.refinementStrategy);

      const auto numFaultElements = meshReader->getFault().size();
      const auto numSubTriangles = faultRefiner->getNumSubTriangles();

      logInfo() << "Initializing Fault output."
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
          constexpr size_t NumVertices{4};
          std::array<const double*, NumVertices> elementVerticesCoords{};
          for (size_t vertexIdx = 0; vertexIdx < NumVertices; ++vertexIdx) {
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
    } else {
      const auto order = elementwiseParams.vtkorder;

      const auto numFaultElements = meshReader->getFault().size();

      // get arrays of elements and vertices from the meshReader
      const auto& faultInfo = meshReader->getFault();
      const auto& elementsInfo = meshReader->getElements();
      const auto& verticesInfo = meshReader->getVertices();

      std::size_t faceCount = 0;
      for (size_t faceIdx = 0; faceIdx < numFaultElements; ++faceIdx) {

        // get a global element ID for the current fault face
        const auto& fault = faultInfo[faceIdx];
        const auto elementIdx = fault.element;

        if (elementIdx >= 0) {
          ++faceCount;
        }
      }

      outputData->receiverPoints.resize(faceCount * seissol::init::vtk2d::Shape[order][1]);
      std::size_t faceOffset = 0;

      // iterate through each fault side
      for (size_t faceIdx = 0; faceIdx < numFaultElements; ++faceIdx) {

        // get a global element ID for the current fault face
        const auto& fault = faultInfo[faceIdx];
        const auto elementIdx = fault.element;

        if (elementIdx >= 0) {
          const auto& element = elementsInfo[elementIdx];

          // store coords of vertices of the current ELEMENT
          constexpr size_t NumVertices{4};
          std::array<const double*, NumVertices> vertices{};
          for (size_t vertexIdx = 0; vertexIdx < NumVertices; ++vertexIdx) {
            auto globalVertexIdx = element.vertices[vertexIdx];
            vertices[vertexIdx] = verticesInfo[globalVertexIdx].coords;
          }

          const auto faceSideIdx = fault.side;

          // init reference coordinates of the fault face
          ExtTriangle referenceTriangle = getReferenceTriangle(faceSideIdx);

          // init global coordinates of the fault face
          ExtTriangle globalFace = getGlobalTriangle(faceSideIdx, element, verticesInfo);

          for (std::size_t i = 0; i < seissol::init::vtk2d::Shape[order][1]; ++i) {
            auto& receiverPoint =
                outputData->receiverPoints[faceOffset * seissol::init::vtk2d::Shape[order][1] + i];
            real nullpoint[2] = {0, 0};
            const real* prepoint =
                i > 0 ? (seissol::init::vtk2d::Values[order] + (i - 1) * 2) : nullpoint;
            double point[2] = {prepoint[0], prepoint[1]};
            transformations::chiTau2XiEtaZeta(faceSideIdx, point, receiverPoint.reference.coords);
            transformations::tetrahedronReferenceToGlobal(vertices[0],
                                                          vertices[1],
                                                          vertices[2],
                                                          vertices[3],
                                                          receiverPoint.reference.coords,
                                                          receiverPoint.global.coords);
            receiverPoint.globalTriangle = globalFace;
            receiverPoint.isInside = true;
            receiverPoint.faultFaceIndex = faceIdx;
            receiverPoint.localFaceSideId = faceSideIdx;
            receiverPoint.elementIndex = element.localId;
            receiverPoint.globalReceiverIndex =
                faceOffset * seissol::init::vtk2d::Shape[order][1] + i;
            receiverPoint.faultTag = fault.tag;
          }

          ++faceOffset;
        }
      }
    }
  }

  inline const static size_t MaxAllowedCacheLevel = 1;

  private:
  seissol::initializer::parameters::ElementwiseFaultParameters elementwiseParams;
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_ELEMENTWISEBUILDER_H_
