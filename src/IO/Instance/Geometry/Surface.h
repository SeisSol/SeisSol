// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_SURFACE_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_SURFACE_H_

#include "Geometry/MeshReader.h"
#include "IO/Instance/Geometry/Geometry.h"

#include <vector>
namespace seissol::io::instance::geometry {

class SurfaceWriter : public GeometryWriter {
  public:
  SurfaceWriter(const std::string& name,
                const ::seissol::geometry::MeshReader& mesh,
                int targetDegree);

  template <typename T>
  void addFaceNodalOutput(const std::string& name,
                          const std::vector<std::size_t>& dimensions,
                          const T* data,
                          std::size_t quantity,
                          std::size_t simulation,
                          bool isConst) {
    std::size_t order(this->order);
    std::vector<real*> locations(indices.size());

    constexpr auto FaceDisplacementPadded =
        tensor::faceDisplacement::Size / tensor::faceDisplacement::Shape[1];
    for (std::size_t i = 0; i < indices.size(); ++i) {
      locations[i] =
          data[meshIndices[i].first][meshIndices[i].second] + FaceDisplacementPadded * quantity;
    }
    addGeometryOutput(
        name, std::vector<std::size_t>(), isConst, [=](real* target, std::size_t index) {
          kernel::projectNodalToVtkFace vtkproj;
          vtkproj.pn = locations[index];
          vtkproj.MV2nTo2m = nodal::init::MV2nTo2m::Values;
          vtkproj.xf(order) = target;
          vtkproj.collff(ConvergenceOrder, order) =
              init::collff::Values[ConvergenceOrder + (ConvergenceOrder + 1) * order];
          vtkproj.execute(order);
        });
  }

  private:
  std::vector<std::size_t> indices;
  std::vector<std::pair<std::size_t, int>> meshIndices;
};

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_SURFACE_H_
