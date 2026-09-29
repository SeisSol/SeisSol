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
#include "Geometry/FaceTransform.h"
#include "Geometry/MeshReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BoundarySetup.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <utils/logger.h>

namespace seissol::initializer {

void initializeBoundaryMappings(const seissol::geometry::MeshReader& meshReader,
                                const std::optional<DirichletCondition>& dirichletCondition,
                                LTS::Storage& ltsStorage) {
  for (auto& layer : ltsStorage.leaves(Ghost)) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* boundary = layer.var<LTS::BoundaryMapping>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      const auto meshId = secondaryInformation[cell].meshId;
      for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
        if (!boundaryProperties(cellInformation[cell].faceTypes[side]).requiresFaceData) {
          continue;
        }

        const auto face =
            seissol::geometry::AffineFaceTransform::fromMeshCell(meshId, side, meshReader);

        // Compute nodal points in global coordinates for each side.
        real nodesReferenceData[nodal::tensor::nodes2D::Size];
        std::copy_n(nodal::init::nodes2D::Values, nodal::tensor::nodes2D::Size, nodesReferenceData);
        auto nodesReference = nodal::init::nodes2D::view::create(nodesReferenceData);
        auto* nodes = boundary[cell][side].nodes;
        assert(nodes != nullptr);
        auto offset = 0;
        for (std::size_t i = 0; i < nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension];
             ++i) {
          // Compute the global coordinates for the nodal points.
          const auto xyz = face.refToSpace(seissol::geometry::FaceTransform::FaceVectorT(
              nodesReference(i, 0), nodesReference(i, 1)));
          for (std::size_t d = 0; d < Cell::Dim; ++d) {
            nodes[offset++] = xyz(d);
          }
        }

        // Compute map that rotates to normal aligned coordinate system.
        real* matTData = boundary[cell][side].dataT;
        real* matTinvData = boundary[cell][side].dataTinv;
        assert(matTData != nullptr);
        assert(matTinvData != nullptr);
        auto matT = init::T::view::create(matTData);
        auto matTinv = init::Tinv::view::create(matTinvData);

        const auto basis = face.faceAlignedBasis();
        seissol::model::getFaceRotationMatrix(
            basis[0].normalized(), basis[1].normalized(), basis[2].normalized(), matT, matTinv);

        // Evaluate easi boundary condition matrices if needed
        real* dirichletMap = boundary[cell][side].dirichletMap;
        real* dirichletOffset = boundary[cell][side].dirichletOffset;
        assert(dirichletMap != nullptr);
        assert(dirichletOffset != nullptr);
        if (cellInformation[cell].faceTypes[side] == FaceType::Dirichlet) {
          if (dirichletCondition.has_value()) {
            const auto faceBarycenter = face.center();

            real globalMapData[tensor::dirichletMapGlobal::size()];
            real globalConstantData[tensor::dirichletOffsetGlobal::size()];
            const auto frame =
                dirichletCondition->query(faceBarycenter.data(), globalMapData, globalConstantData);

            if (frame == BoundaryFrame::FaceAligned) {
              std::copy_n(globalMapData, tensor::dirichletMap::size(), dirichletMap);
              std::copy_n(globalConstantData, tensor::dirichletOffset::size(), dirichletOffset);
            } else {
              kernel::rotateBoundaryCondition rotateKrnl;
              rotateKrnl.dirichletMapGlobal = globalMapData;
              rotateKrnl.dirichletOffsetGlobal = globalConstantData;
              rotateKrnl.dirichletMap = dirichletMap;
              rotateKrnl.dirichletOffset = dirichletOffset;
              rotateKrnl.T = matTData;
              rotateKrnl.Tinv = matTinvData;
              rotateKrnl.execute();
            }
          } else {
            logError() << "Dirichlet face found, but no boundary condition definition given.";
          }
        } else {
          // Boundary should not be evaluated
          std::fill_n(dirichletMap,
                      seissol::tensor::dirichletMap::size(),
                      std::numeric_limits<real>::signaling_NaN());
          std::fill_n(dirichletOffset,
                      seissol::tensor::dirichletOffset::size(),
                      std::numeric_limits<real>::signaling_NaN());
        }
      }
    }
  }
}

} // namespace seissol::initializer
