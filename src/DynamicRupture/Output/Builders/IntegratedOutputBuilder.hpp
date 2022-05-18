#ifndef SEISSOL_DR_INTEGRATED_OUTPUT_BUILDER_HPP
#define SEISSOL_DR_INTEGRATED_OUTPUT_BUILDER_HPP

#include "Geometry/MeshReader.h"
#include "Initializer/tree/Lut.hpp"
#include "Initializer/LTS.h"

namespace seissol::dr::output {
class IntegratedOutputBuilder {
  public:
  void setMeshReader(const MeshReader* reader) { meshReader = reader; }
  void setLtsData(seissol::initializers::LTSTree* userWpTree,
                  seissol::initializers::LTS* userWpDescr,
                  seissol::initializers::Lut* userWpLut) {
    wpTree = userWpTree;
    wpDescr = userWpDescr;
    wpLut = userWpLut;
  }

  void build(IntegratedOutputData* geoOutputData) {
    outputData = geoOutputData;
    computeSurfaceAreas();
    saveMaterialData();
  }

  protected:
  void computeSurfaceAreas() {
    const auto numFaultElements = meshReader->getFault().size();

    // get arrays of elements and vertices from the meshReader
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    // iterate through each fault side
    for (size_t faceIndex = 0; faceIndex < numFaultElements; ++faceIndex) {

      // get a global element ID for the current fault face
      auto elementIndex = faultInfo[faceIndex].element;
      if (elementIndex >= 0) {
        const auto& element = elementsInfo[elementIndex];
        auto localFaceSideId = faultInfo[faceIndex].side;

        // find global coordinates of the fault face
        ExtTriangle face = getGlobalTriangle(localFaceSideId, element, verticesInfo);

        // compute area and append the vector
        outputData->surfaceAreas.push_back(computeTriangleArea(face));
      } else {
        outputData->surfaceAreas.push_back(0.0);
      }
    }
  }

  void saveMaterialData() {
    const auto numFaultElements = meshReader->getFault().size();

    // get arrays of elements and vertices from the meshReader
    const auto& faultInfos = meshReader->getFault();

    // iterate through each fault side
    for (size_t faceIndex = 0; faceIndex < numFaultElements; ++faceIndex) {
      auto faultInfo = faultInfos[faceIndex];
      auto elementIndex = faultInfos[faceIndex].element;
      if (elementIndex >= 0) {
        auto& local = wpLut->lookup(wpDescr->material, faultInfo.element).local;

#if defined USE_ANISOTROPIC
        double muBar = (local.c44 + local.c55 + local.c66) / 3.0;
        auto lambda = (local.c11 + local.c22 + local.c33) / 3.0 - 2.0 * muBar;
#else
        auto lambda = local.lambda;
#endif
        outputData->lambda.push_back(lambda);
      } else {
        outputData->lambda.push_back(0.0);
      }
    }
  }

  private:
  const MeshReader* meshReader{};
  IntegratedOutputData* outputData{};

  seissol::initializers::LTS* wpDescr{nullptr};
  seissol::initializers::LTSTree* wpTree{nullptr};
  seissol::initializers::Lut* wpLut{nullptr};
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_INTEGRATED_OUTPUT_BUILDER_HPP
