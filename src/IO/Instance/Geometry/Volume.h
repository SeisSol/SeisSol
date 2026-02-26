// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_VOLUME_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_VOLUME_H_

#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Instance/Geometry/Points.h"
#include "IO/Instance/Geometry/Refinement.h"

#include <vector>
namespace seissol::io::instance::geometry {

enum class RefinementMode : int32_t { Refine1 = 1, Refine4 = 4, Refine8 = 8, Refine32 = 32 };

class VolumeWriter : public GeometryWriter {
  public:
  VolumeWriter(const std::string& name,
               const ::seissol::geometry::MeshReader& mesh,
               const WriterConfig& config,
               RefinementMode refinement)
      : GeometryWriter(name,
                       mesh.getElements().size(),
                       Shape::Tetrahedron,
                       config,
                       static_cast<int32_t>(refinement)) {

    const auto dataOrder = config.order;
    const auto trueOrder = std::max(config.order, 1U);
    const auto trueBase = pointsTetrahedron(trueOrder);
    const auto dataBase = pointsTetrahedron(dataOrder);

    auto truePoints = std::vector<std::vector<std::array<double, 3>>>{trueBase};
    auto dataPoints = std::vector<std::vector<std::array<double, 3>>>{dataBase};

    if (refinement == RefinementMode::Refine4) {
      truePoints = applySubdivide(truePoints, TetrahedronRefine4);
      dataPoints = applySubdivide(dataPoints, TetrahedronRefine4);
    }
    if (refinement == RefinementMode::Refine8) {
      truePoints = applySubdivide(truePoints, TetrahedronRefine8);
      dataPoints = applySubdivide(dataPoints, TetrahedronRefine8);
    }
    if (refinement == RefinementMode::Refine32) {
      truePoints = applySubdivide(truePoints, TetrahedronRefine4);
      dataPoints = applySubdivide(dataPoints, TetrahedronRefine4);
      truePoints = applySubdivide(truePoints, TetrahedronRefine8);
      dataPoints = applySubdivide(dataPoints, TetrahedronRefine8);
    }

    const auto& elements = mesh.getElements();
    const auto& vertices = mesh.getVertices();
    addPointProjector([=](double* target, std::size_t cell, std::size_t subcell) {
      const auto& elemVertices = elements[cell].vertices;

      for (std::size_t i = 0; i < truePoints[subcell].size(); ++i) {
        seissol::transformations::tetrahedronReferenceToGlobal(vertices[elemVertices[0]].coords,
                                                               vertices[elemVertices[1]].coords,
                                                               vertices[elemVertices[2]].coords,
                                                               vertices[elemVertices[3]].coords,
                                                               truePoints[subcell][i].data(),
                                                               &target[i * Cell::Dim]);
      }
    });

    for (const auto& tetrahedron : dataPoints) {
      auto& coll = proj.emplace_back();
      std::size_t idx = 0;
      for (const auto& point : tetrahedron) {
        const auto data = basisFunction::SampledBasisFunctions<real>(
                              ConvergenceOrder, point[0], point[1], point[2])
                              .m_data;
        std::copy(data.begin(), data.end(), coll.begin() + idx);
        idx += data.size();
      }
    }
  }

  template <typename T>
  void addVolumeBasisOutput(const std::string& name,
                            const std::vector<const T*>& data,
                            std::size_t simulation,
                            bool isConst) {
    std::size_t order(this->config.order);
    const auto proj = this->proj;

    addGeometryOutput(name,
                      std::vector<std::size_t>(),
                      isConst,
                      [=](real* target, std::size_t cell, std::size_t subcell) {
                        if (data[cell] != nullptr) {
                          kernel::projectBasisToVtkVolume vtkproj;
                          vtkproj.qb = data[cell];
                          vtkproj.xv(order) = target;
                          vtkproj.collvv(ConvergenceOrder, order) = proj[subcell].data();
                          vtkproj.execute(order);
                        } else {
                          std::memset(target, 0, sizeof(real)); // TODO
                        }
                      });
  }

  private:
  std::vector<std::size_t> indices;
  std::vector<
      memory::AlignedArray<real, static_cast<std::size_t>(tensor::Q::Size) * tensor::Q::Size>>
      proj;
};

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_VOLUME_H_
