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
#include "Geometry/FaceTransform.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "ReceiverBasedOutputBuilder.h"
#include "Solver/MultipleSimulations.h"

namespace seissol::dr::output {
class ElementWiseBuilder : public ReceiverBasedOutputBuilder {
  public:
  ~ElementWiseBuilder() override = default;
  void setParams(const seissol::initializer::parameters::ElementwiseFaultParameters& params) {
    elementwiseParams_ = params;
  }
  void build(std::shared_ptr<ReceiverOutputData> elementwiseOutputData) {
    outputData_ = std::move(elementwiseOutputData);
    initReceiverLocations();
    assignNearestGaussianPoints(outputData_->receivers);
    assignNearestInternalGaussianPoints();
    assignFusedIndices();
    assignFaultTags();
    initTimeCaching();
    // initTopology establishes the face/point hierarchy all following steps index into, and fixes
    // the receiver numbering; everything below has to run after it
    initTopology();
    initOutputVariables(elementwiseParams_.outputMask);
    initBasisFunctions();
    initDeviceCollectors(true);
    initFaultDirections();
    initRotationMatrices();
    initJacobian2dMatrices();
    outputData_->isActive = true;
  }

  protected:
  void initTimeCaching() override {
    outputData_->maxCacheLevel = ElementWiseBuilder::MaxAllowedCacheLevel;
    outputData_->currentCacheLevel = 0;
  }

  void initReceiverLocations() {
    auto faultRefiner = refiner::get(elementwiseParams_.refinementStrategy);

    const auto numSubTriangles = faultRefiner->getNumSubTriangles();
    const auto order = static_cast<std::uint32_t>(std::max(elementwiseParams_.vtkorder, 0));

    const auto numFaultElements = meshReader_->getFault().size();

    logInfo() << "Initializing Fault output."
              << "Number of sub-triangles:" << numSubTriangles << "Output order:" << order
              << "Simulation count:" << multisim::NumSimulations;

    // get the array of fault faces from the meshReader
    const auto& faultInfo = meshReader_->getFault();

    // iterate through each fault side
    for (size_t faceIdx = 0; faceIdx < numFaultElements; ++faceIdx) {
      const auto& fault = faultInfo[faceIdx];
      const auto elementIdx = fault.element;

      if (elementIdx.hasValue()) {
        const auto faceSideIdx = fault.side;

        // init reference coordinates of the fault face
        const ExtTriangle referenceTriangle = getReferenceTriangle(faceSideIdx);

        // init global coordinates of the fault face
        const ExtTriangle globalFace =
            toExtTriangle(seissol::geometry::AffineFaceTransform::fromMeshCell(
                elementIdx.value(), faceSideIdx, *meshReader_));

        faultRefiner->refineAndAccumulate({elementwiseParams_.refinement,
                                           faceIdx,
                                           faceSideIdx,
                                           elementIdx.value(),
                                           &fault,
                                           order,
                                           multisim::NumSimulations},
                                          std::make_pair(globalFace, referenceTriangle));
      }
    }

    // retrieve all receivers from a fault face refiner
    outputData_->receivers = faultRefiner->moveAllReceivers();
    faultRefiner.reset(nullptr);

    // the refiner divides the plane face through the vertices, at the same reference coordinates
    // as the face of the mesh
    placeOnCurvedFaces(outputData_->receivers);
  }

  inline const static size_t MaxAllowedCacheLevel = 1;

  private:
  seissol::initializer::parameters::ElementwiseFaultParameters elementwiseParams_;
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_ELEMENTWISEBUILDER_H_
