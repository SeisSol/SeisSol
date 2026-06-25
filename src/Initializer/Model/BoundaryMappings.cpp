// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "BoundaryMappings.h"

#include "Common/Constants.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <vector>

namespace seissol::initializer {

void initializeBoundaryMappings(const seissol::geometry::MeshReader& meshReader,
                                const std::optional<EasiBoundary>& easiBoundary,
                                LTS::Storage& ltsStorage) {
  const std::vector<Element>& elements = meshReader.getElements();
  const std::vector<Vertex>& vertices = meshReader.getVertices();

  for (auto& layer : ltsStorage.leaves(Ghost)) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* boundary = layer.var<LTS::BoundaryMapping>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      const auto& element = elements[secondaryInformation[cell].meshId];
      const double* coords[Cell::NumVertices];
      for (std::size_t v = 0; v < Cell::NumVertices; ++v) {
        coords[v] = vertices[element.vertices[v]].coords;
      }
      for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
        if (cellInformation[cell].faceTypes[side] != FaceType::FreeSurfaceGravity &&
            cellInformation[cell].faceTypes[side] != FaceType::Dirichlet &&
            cellInformation[cell].faceTypes[side] != FaceType::Analytical) {
          continue;
        }
        // Compute nodal points in global coordinates for each side.
        real nodesReferenceData[nodal::tensor::nodes2D::Size];
        std::copy_n(nodal::init::nodes2D::Values, nodal::tensor::nodes2D::Size, nodesReferenceData);
        auto nodesReference = nodal::init::nodes2D::view::create(nodesReferenceData);
        auto* nodes = boundary[cell][side].nodes;
        assert(nodes != nullptr);
        auto offset = 0;
        for (std::size_t i = 0; i < nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension];
             ++i) {
          double nodeReference[2];
          nodeReference[0] = nodesReference(i, 0);
          nodeReference[1] = nodesReference(i, 1);
          // Compute the global coordinates for the nodal points.
          double xiEtaZeta[3];
          double xyz[3];
          seissol::transformations::chiTau2XiEtaZeta(side, nodeReference, xiEtaZeta);
          seissol::transformations::tetrahedronReferenceToGlobal(
              coords[0], coords[1], coords[2], coords[3], xiEtaZeta, xyz);
          nodes[offset++] = xyz[0];
          nodes[offset++] = xyz[1];
          nodes[offset++] = xyz[2];
        }

        // Compute map that rotates to normal aligned coordinate system.
        real* matTData = boundary[cell][side].dataT;
        real* matTinvData = boundary[cell][side].dataTinv;
        assert(matTData != nullptr);
        assert(matTinvData != nullptr);
        auto matT = init::T::view::create(matTData);
        auto matTinv = init::Tinv::view::create(matTinvData);

        VrtxCoords normal;
        VrtxCoords tangent1;
        VrtxCoords tangent2;
        MeshTools::normalAndTangents(element, side, vertices, normal, tangent1, tangent2);
        MeshTools::normalize(normal, normal);
        MeshTools::normalize(tangent1, tangent1);
        MeshTools::normalize(tangent2, tangent2);
        seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, matT, matTinv);

        // Evaluate easi boundary condition matrices if needed
        real* easiBoundaryMap = boundary[cell][side].easiBoundaryMap;
        real* easiBoundaryConstant = boundary[cell][side].easiBoundaryConstant;
        assert(easiBoundaryMap != nullptr);
        assert(easiBoundaryConstant != nullptr);
        if (cellInformation[cell].faceTypes[side] == FaceType::Dirichlet) {
          if (easiBoundary.has_value()) {
            easiBoundary->query(nodes, easiBoundaryMap, easiBoundaryConstant);
          } else {
            logError() << "Dirichlet face found, but no boundary condition definition given.";
          }
        } else {
          // Boundary should not be evaluated
          std::fill_n(easiBoundaryMap,
                      seissol::tensor::easiBoundaryMap::size(),
                      std::numeric_limits<real>::signaling_NaN());
          std::fill_n(easiBoundaryConstant,
                      seissol::tensor::easiBoundaryConstant::size(),
                      std::numeric_limits<real>::signaling_NaN());
        }
      }
    }
  }
}

} // namespace seissol::initializer
