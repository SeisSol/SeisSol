#ifndef SEISSOL_GEOMETRY_BUILDER_HPP
#define SEISSOL_GEOMETRY_BUILDER_HPP

#include "Geometry/MeshReader.h"

namespace seissol::dr::output {
class GeometryBuilder {
  public:
  void setMeshReader(const MeshReader* reader) { meshReader = reader; }

  void build(GeoOutputData &outputData) {
    const auto numFaultElements = meshReader->getFault().size();

    // get arrays of elements and vertices from the meshReader
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    // iterate through each fault side
    for (size_t faceIndex = 0; faceIndex < numFaultElements; ++faceIndex) {

      // get a global element ID for the current fault face
      auto elementIndex = faultInfo[faceIndex].element;
      const auto& element = elementsInfo[elementIndex];
      auto localFaceSideId = faultInfo[faceIndex].side;

      // find global coordinates of the fault face
      ExtTriangle face = getGlobalTriangle(localFaceSideId, element, verticesInfo);

      // compute area and append the vector
      outputData.surfaceAreas.push_back(computeTriangleArea(face));
    }
  }

  private:
  const MeshReader* meshReader{};
};
} // namespace seissol::dr::output

#endif // SEISSOL_GEOMETRY_BUILDER_HPP
