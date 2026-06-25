// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "CellLocalMatrices.h"

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
#include "Initializer/MemoryManager.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Model/CommonDatastructures.h"
#include "Numerical/Transformation.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace seissol::initializer {

namespace {

void setStarMatrix(const real* matAT,
                   const real* matBT,
                   const real* matCT,
                   const double grad[3],
                   real* starMatrix) {
  for (std::size_t idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
    starMatrix[idx] = grad[0] * matAT[idx];
  }

  for (std::size_t idx = 0; idx < seissol::tensor::star::size(1); ++idx) {
    starMatrix[idx] += grad[1] * matBT[idx];
  }

  for (std::size_t idx = 0; idx < seissol::tensor::star::size(2); ++idx) {
    starMatrix[idx] += grad[2] * matCT[idx];
  }
}

} // namespace

void initializeCellLocalMatrices(const seissol::geometry::MeshReader& meshReader,
                                 LTS::Storage& ltsStorage,
                                 const ClusterLayout& clusterLayout,
                                 const parameters::ModelParameters& modelParameters) {
  const std::vector<Element>& elements = meshReader.getElements();
  const std::vector<Vertex>& vertices = meshReader.getVertices();

  static_assert(seissol::tensor::AplusT::Shape[0] == seissol::tensor::AminusT::Shape[0],
                "Shape mismatch for flux matrices");
  static_assert(seissol::tensor::AplusT::Shape[1] == seissol::tensor::AminusT::Shape[1],
                "Shape mismatch for flux matrices");

  assert(LayerMask(Ghost) == ltsStorage.info<LTS::Material>().mask);
  assert(LayerMask(Ghost) == ltsStorage.info<LTS::LocalIntegration>().mask);
  assert(LayerMask(Ghost) == ltsStorage.info<LTS::NeighboringIntegration>().mask);

  for (auto& layer : ltsStorage.leaves(Ghost)) {
    auto* material = layer.var<LTS::Material>();
    auto* materialData = layer.var<LTS::MaterialData>();
    auto* localIntegration = layer.var<LTS::LocalIntegration>();
    auto* neighboringIntegration = layer.var<LTS::NeighboringIntegration>();
    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();

#pragma omp parallel
    {
      real matATData[tensor::star::size(0)];
      real matATtildeData[tensor::star::size(0)];
      real matBTData[tensor::star::size(1)];
      real matCTData[tensor::star::size(2)];
      auto matAT = init::star::view<0>::create(matATData);
      // matAT with elastic parameters in local coordinate system, used for flux kernel
      auto matATtilde = init::star::view<0>::create(matATtildeData);
      auto matBT = init::star::view<0>::create(matBTData);
      auto matCT = init::star::view<0>::create(matCTData);

      real matTData[seissol::tensor::T::size()];
      real matTinvData[seissol::tensor::Tinv::size()];
      auto matT = init::T::view::create(matTData);
      auto matTinv = init::Tinv::view::create(matTinvData);

      real qGodLocalData[tensor::QgodLocal::size()];
      real qGodNeighborData[tensor::QgodNeighbor::size()];
      auto qGodLocal = init::QgodLocal::view::create(qGodLocalData);
      auto qGodNeighbor = init::QgodNeighbor::view::create(qGodNeighborData);

      real rusanovPlusNull[tensor::QcorrLocal::size()]{};
      real rusanovMinusNull[tensor::QcorrNeighbor::size()]{};

#pragma omp for schedule(static)
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        const auto clusterId = secondaryInformation[cell].clusterId;
        const auto timeStepWidth = clusterLayout.timestepRate(clusterId);
        const auto meshId = secondaryInformation[cell].meshId;

        // NOLINTNEXTLINE
        auto& materialLocal = materialData[cell];

        double x[Cell::NumVertices];
        double y[Cell::NumVertices];
        double z[Cell::NumVertices];
        double gradXi[3];
        double gradEta[3];
        double gradZeta[3];

        // Iterate over all 4 vertices of the tetrahedron
        for (std::size_t vertex = 0; vertex < Cell::NumVertices; ++vertex) {
          const VrtxCoords& coords = vertices[elements[meshId].vertices[vertex]].coords;
          x[vertex] = coords[0];
          y[vertex] = coords[1];
          z[vertex] = coords[2];
        }

        seissol::transformations::tetrahedronGlobalToReferenceJacobian(
            x, y, z, gradXi, gradEta, gradZeta);

        seissol::model::getTransposedCoefficientMatrix(materialLocal, 0, matAT);
        seissol::model::getTransposedCoefficientMatrix(materialLocal, 1, matBT);
        seissol::model::getTransposedCoefficientMatrix(materialLocal, 2, matCT);

        setStarMatrix(
            matATData, matBTData, matCTData, gradXi, localIntegration[cell].starMatrices[0]);
        setStarMatrix(
            matATData, matBTData, matCTData, gradEta, localIntegration[cell].starMatrices[1]);
        setStarMatrix(
            matATData, matBTData, matCTData, gradZeta, localIntegration[cell].starMatrices[2]);

        const double volume = MeshTools::volume(elements[meshId], vertices);

        for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
          VrtxCoords normal;
          VrtxCoords tangent1;
          VrtxCoords tangent2;
          MeshTools::normalAndTangents(
              elements[meshId], side, vertices, normal, tangent1, tangent2);
          const double surface = MeshTools::surface(normal);
          MeshTools::normalize(normal, normal);
          MeshTools::normalize(tangent1, tangent1);
          MeshTools::normalize(tangent2, tangent2);

          // Defines a rotation matrix for computing material properties in face-local coordinates
          // for anisotropy. It has no effect for isotropic materials.
          std::array<double, 36> nLocalData{};
          seissol::model::getBondMatrix(normal, tangent1, tangent2, nLocalData);
          seissol::model::getTransposedGodunovState(
              seissol::model::getRotatedMaterialCoefficients(nLocalData, materialLocal),
              seissol::model::getRotatedMaterialCoefficients(
                  nLocalData, *dynamic_cast<model::MaterialT*>(material[cell].neighbor[side])),
              cellInformation[cell].faceTypes[side],
              qGodLocal,
              qGodNeighbor);
          seissol::model::getTransposedCoefficientMatrix(
              seissol::model::getRotatedMaterialCoefficients(nLocalData, materialLocal),
              0,
              matATtilde);

          // Calculate transposed T and Tinv instead
          seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, matT, matTinv);

          // Scale with |S_side|/|J| and multiply with -1 as the flux matrices
          // must be subtracted.
          const double fluxScale = -2.0 * surface / (6.0 * volume);

          const auto isSpecialBC = [&](int side) {
            const auto hasDRFace = [](const CellLocalInformation& ci) {
              bool hasAtLeastOneDRFace = false;
              for (size_t i = 0; i < Cell::NumFaces; ++i) {
                if (ci.faceTypes[i] == FaceType::DynamicRupture) {
                  hasAtLeastOneDRFace = true;
                }
              }
              return hasAtLeastOneDRFace;
            };
            const bool thisCellHasAtLeastOneDRFace = hasDRFace(cellInformation[cell]);
            const auto& neighborID = secondaryInformation[cell].faceNeighbors[side];
            const bool neighborBehindSideHasAtLeastOneDRFace =
                neighborID != StoragePosition::NullPosition &&
                hasDRFace(ltsStorage.lookup<LTS::CellInformation>(neighborID));
            const bool adjacentDRFaceExists =
                thisCellHasAtLeastOneDRFace || neighborBehindSideHasAtLeastOneDRFace;
            return (cellInformation[cell].faceTypes[side] == FaceType::Regular) &&
                   adjacentDRFaceExists;
          };

          const auto wavespeedLocal = materialLocal.getMaxWaveSpeed();
          const auto wavespeedNeighbor = material[cell].neighbor[side]->getMaxWaveSpeed();
          const auto wavespeed = std::max(wavespeedLocal, wavespeedNeighbor);

          real centralFluxData[tensor::QgodLocal::size()]{};
          real rusanovPlusData[tensor::QcorrLocal::size()]{};
          real rusanovMinusData[tensor::QcorrNeighbor::size()]{};
          auto centralFluxView = init::QgodLocal::view::create(centralFluxData);
          auto rusanovPlusView = init::QcorrLocal::view::create(rusanovPlusData);
          auto rusanovMinusView = init::QcorrNeighbor::view::create(rusanovMinusData);
          for (size_t i = 0; i < std::min(tensor::QgodLocal::Shape[0], tensor::QgodLocal::Shape[1]);
               i++) {
            centralFluxView(i, i) = 0.5;
            rusanovPlusView(i, i) = wavespeed * 0.5;
            rusanovMinusView(i, i) = -wavespeed * 0.5;
          }

          // check if we're on a face that has an adjacent cell with DR face
          const auto fluxDefault =
              isSpecialBC(side) ? modelParameters.fluxNearFault : modelParameters.flux;

          // exclude boundary conditions
          static const std::vector<FaceType> GodunovBoundaryConditions = {
              FaceType::FreeSurface,
              FaceType::FreeSurfaceGravity,
              FaceType::Analytical,
              FaceType::Outflow};

          const auto enforceGodunovBc = std::any_of(
              GodunovBoundaryConditions.begin(),
              GodunovBoundaryConditions.end(),
              [&](auto condition) { return condition == cellInformation[cell].faceTypes[side]; });

          const auto enforceGodunovEa = isAtElasticAcousticInterface(material[cell], side);

          const auto enforceGodunov = enforceGodunovBc || enforceGodunovEa;

          const auto flux = enforceGodunov ? parameters::NumericalFlux::Godunov : fluxDefault;

          kernel::computeFluxSolverLocal localKrnl;
          localKrnl.fluxScale = fluxScale;
          localKrnl.AplusT = localIntegration[cell].nApNm1[side];
          if (cellInformation[cell].faceTypes[side] == FaceType::DynamicRupture) {
            localKrnl.fluxScale = 0;
          }
          if (flux == parameters::NumericalFlux::Rusanov) {
            localKrnl.QgodLocal = centralFluxData;
            localKrnl.QcorrLocal = rusanovPlusData;
          } else {
            localKrnl.QgodLocal = qGodLocalData;
            localKrnl.QcorrLocal = rusanovPlusNull;
          }
          localKrnl.T = matTData;
          localKrnl.Tinv = matTinvData;
          localKrnl.star(0) = matATtildeData;
          localKrnl.execute();

          kernel::computeFluxSolverNeighbor neighKrnl;
          neighKrnl.fluxScale = fluxScale;
          neighKrnl.AminusT = neighboringIntegration[cell].nAmNm1[side];
          if (flux == parameters::NumericalFlux::Rusanov) {
            neighKrnl.QgodNeighbor = centralFluxData;
            neighKrnl.QcorrNeighbor = rusanovMinusData;
          } else {
            neighKrnl.QgodNeighbor = qGodNeighborData;
            neighKrnl.QcorrNeighbor = rusanovMinusNull;
          }
          neighKrnl.T = matTData;
          neighKrnl.Tinv = matTinvData;
          neighKrnl.star(0) = matATtildeData;
          if (cellInformation[cell].faceTypes[side] == FaceType::Dirichlet ||
              cellInformation[cell].faceTypes[side] == FaceType::FreeSurfaceGravity) {
            // already rotated
            neighKrnl.Tinv = init::identityT::Values;
          }
          neighKrnl.execute();
        }

        seissol::model::initializeSpecificLocalData(
            materialLocal, timeStepWidth, &localIntegration[cell].specific);

        seissol::model::initializeSpecificNeighborData(materialLocal,
                                                       &neighboringIntegration[cell].specific);
      }
    }
  }
}

} // namespace seissol::initializer
