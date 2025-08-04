// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_VOLUME_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_VOLUME_H_

#include <IO/Instance/Geometry/Geometry.h>
#include <vector>
namespace seissol::io::instance::geometry {

class VolumeWriter : public GeometryWriter {
  public:
  VolumeWriter(const std::string& name,
               const ::seissol::geometry::MeshReader& mesh,
               int targetDegree);

  template <typename T>
  void addVolumeBasisOutput(const std::string& name,
                            const std::vector<std::size_t>& dimensions,
                            const T* data,
                            std::size_t quantity,
                            std::size_t simulation,
                            bool isConst) {
    std::size_t order(this->order);
    std::vector<real*> locations(indices.size());

    constexpr auto QDofSizePadded = tensor::Q::Size / tensor::Q::Shape[1];
    for (std::size_t i = 0; i < indices.size(); ++i) {
      locations[i] = data[indices[i]] + QDofSizePadded * quantity;
    }
    addGeometryOutput(
        name, std::vector<std::size_t>(), isConst, [=](real* target, std::size_t index) {
          kernel::projectBasisToVtkVolume vtkproj;
          vtkproj.qb = locations[index];
          vtkproj.xv(order) = target;
          vtkproj.collvv(ConvergenceOrder, order) =
              init::collvv::Values[ConvergenceOrder + (ConvergenceOrder + 1) * order];
          vtkproj.execute(order);
        });
  }

  private:
  std::vector<std::size_t> indices;
};

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_VOLUME_H_
