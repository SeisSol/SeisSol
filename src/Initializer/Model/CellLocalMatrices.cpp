// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "CellLocalMatrices.h"

#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Config.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BoundaryHelper.h"
#include "Initializer/BoundarySetup.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Model/CommonDatastructures.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::initializer {

namespace {

void setStarMatrix(const real* matAT,
                   const real* matBT,
                   const real* matCT,
                   const std::array<double, Cell::Dim>& grad,
                   real* starMatrix) {
  for (std::size_t idx = 0; idx < seissol::tensor::star<Config>::size(0); ++idx) {
    starMatrix[idx] = grad[0] * matAT[idx];
  }

  for (std::size_t idx = 0; idx < seissol::tensor::star<Config>::size(1); ++idx) {
    starMatrix[idx] += grad[1] * matBT[idx];
  }

  for (std::size_t idx = 0; idx < seissol::tensor::star<Config>::size(2); ++idx) {
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
  constexpr auto Variant = configIdOf<Config>();

  static_assert(seissol::tensor::AplusT<Config>::Shape[0] ==
                    seissol::tensor::AminusT<Config>::Shape[0],
                "Shape mismatch for flux matrices");
  static_assert(seissol::tensor::AplusT<Config>::Shape[1] ==
                    seissol::tensor::AminusT<Config>::Shape[1],
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
    auto* boundaryMapping = layer.var<LTS::BoundaryMapping>();

#pragma omp parallel
    {
      real matATData[tensor::star<Config>::size(0)]{};
      real matATtildeData[tensor::star<Config>::size(0)]{};
      real matBTData[tensor::star<Config>::size(1)]{};
      real matCTData[tensor::star<Config>::size(2)]{};
      auto matAT = init::star<Config>::view<0>::create(matATData);
      // matAT with elastic parameters in local coordinate system, used for flux kernel
      auto matATtilde = init::star<Config>::view<0>::create(matATtildeData);
      auto matBT = init::star<Config>::view<0>::create(matBTData);
      auto matCT = init::star<Config>::view<0>::create(matCTData);

      real matTData[seissol::tensor::T<Config>::size()]{};
      real matTinvData[seissol::tensor::Tinv<Config>::size()]{};
      auto matT = init::T<Config>::view::create(matTData);
      auto matTinv = init::Tinv<Config>::view::create(matTinvData);

      real qGodLocalData[tensor::QgodLocal<Config>::size()]{};
      real qGodNeighborData[tensor::QgodNeighbor<Config>::size()]{};
      auto qGodLocal = init::QgodLocal<Config>::view::create(qGodLocalData);
      auto qGodNeighbor = init::QgodNeighbor<Config>::view::create(qGodNeighborData);

      real rusanovPlusNull[tensor::QcorrLocal<Config>::size()]{};
      real rusanovMinusNull[tensor::QcorrNeighbor<Config>::size()]{};

#pragma omp for schedule(static)
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        const auto clusterId = secondaryInformation[cell].clusterId;
        const auto timeStepWidth = clusterLayout.timestepRate(clusterId);
        const auto meshId = secondaryInformation[cell].meshId;

        // NOLINTNEXTLINE
        auto& materialLocal = materialData[cell];

        std::array<double, Cell::Dim> gradXi{};
        std::array<double, Cell::Dim> gradEta{};
        std::array<double, Cell::Dim> gradZeta{};

        const auto transform = seissol::geometry::AffineTransform::fromMeshCell(meshId, meshReader);

        // IMPORTANT NOTE: we rely on the linearity of the cell transform in this place.
        // hence, you may use an AffineTransform with an arbitrary point here; but nothing more.
        const auto grad = transform.refToSpaceJacobianInverse(
            seissol::geometry::CellTransform::VectorEigenT(Cell::ReferenceBarycenter.data()));

        for (std::size_t i = 0; i < Cell::Dim; ++i) {
          gradXi[i] = grad(0, i);
          gradEta[i] = grad(1, i);
          gradZeta[i] = grad(2, i);
        }

        seissol::model::getTransposedCoefficientMatrix<Config>(materialLocal, 0, matAT);
        seissol::model::getTransposedCoefficientMatrix<Config>(materialLocal, 1, matBT);
        seissol::model::getTransposedCoefficientMatrix<Config>(materialLocal, 2, matCT);

        setStarMatrix(
            matATData, matBTData, matCTData, gradXi, localIntegration[cell].starMatrices[0]);
        setStarMatrix(
            matATData, matBTData, matCTData, gradEta, localIntegration[cell].starMatrices[1]);
        setStarMatrix(
            matATData, matBTData, matCTData, gradZeta, localIntegration[cell].starMatrices[2]);

        const double volume = MeshTools::volume(elements[meshId], vertices);

        for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
          CoordinateT normal{};
          CoordinateT tangent1{};
          CoordinateT tangent2{};
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
          seissol::model::getTransposedCoefficientMatrix<Config>(
              seissol::model::getRotatedMaterialCoefficients(nLocalData, materialLocal),
              0,
              matATtilde);

          // Calculate transposed T and Tinv instead
          seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, matT, matTinv);

          // Scale with |S_side|/|J| and multiply with -1 as the flux matrices
          // must be subtracted.
          const double fluxScale = -2.0 * surface / (6.0 * volume);

          const auto isSpecialBC = [&](std::int8_t side) {
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

          real centralFluxData[tensor::QgodLocal<Config>::size()]{};
          real rusanovPlusData[tensor::QcorrLocal<Config>::size()]{};
          real rusanovMinusData[tensor::QcorrNeighbor<Config>::size()]{};
          auto centralFluxView = init::QgodLocal<Config>::view::create(centralFluxData);
          auto rusanovPlusView = init::QcorrLocal<Config>::view::create(rusanovPlusData);
          auto rusanovMinusView = init::QcorrNeighbor<Config>::view::create(rusanovMinusData);
          for (size_t i = 0; i < std::min(tensor::QgodLocal<Config>::Shape[0],
                                          tensor::QgodLocal<Config>::Shape[1]);
               i++) {
            centralFluxView(i, i) = 0.5;
            rusanovPlusView(i, i) = wavespeed * 0.5;
            rusanovMinusView(i, i) = -wavespeed * 0.5;
          }

          // check if we're on a face that has an adjacent cell with DR face
          const auto fluxDefault =
              isSpecialBC(side) ? modelParameters.fluxNearFault : modelParameters.flux;

          const auto enforceGodunovBc =
              boundaryProperties(cellInformation[cell].faceTypes[side]).enforcesGodunovFlux;

          const auto enforceGodunovEa = isAtElasticAcousticInterface(material[cell], side);

          const auto enforceGodunov = enforceGodunovBc || enforceGodunovEa;

          const auto flux = enforceGodunov ? parameters::NumericalFlux::Godunov : fluxDefault;

          runtime::kernel::computeFluxSolverLocal localKrnl;
          localKrnl.fluxScale = fluxScale;
          localKrnl.AplusT =
              runtime::init::AplusT::view(Variant, localIntegration[cell].nApNm1[side]);
          if (cellInformation[cell].faceTypes[side] == FaceType::DynamicRupture) {
            localKrnl.fluxScale = 0;
          }
          if (flux == parameters::NumericalFlux::Rusanov) {
            localKrnl.QgodLocal = runtime::init::QgodLocal::view(Variant, centralFluxData);
            localKrnl.QcorrLocal = runtime::init::QcorrLocal::view(Variant, rusanovPlusData);
          } else {
            localKrnl.QgodLocal = runtime::init::QgodLocal::view(Variant, qGodLocalData);
            localKrnl.QcorrLocal = runtime::init::QcorrLocal::view(Variant, rusanovPlusNull);
          }
          localKrnl.T = runtime::init::T::view(Variant, matTData);
          localKrnl.Tinv = runtime::init::Tinv::view(Variant, matTinvData);
          localKrnl.star(0) = runtime::init::star::view(Variant, 0, matATtildeData);
          localKrnl.execute(Variant);

          runtime::kernel::computeFluxSolverNeighbor neighKrnl;
          neighKrnl.fluxScale = fluxScale;
          neighKrnl.AminusT =
              runtime::init::AminusT::view(Variant, neighboringIntegration[cell].nAmNm1[side]);
          if (flux == parameters::NumericalFlux::Rusanov) {
            neighKrnl.QgodNeighbor = runtime::init::QgodNeighbor::view(Variant, centralFluxData);
            neighKrnl.QcorrNeighbor = runtime::init::QcorrNeighbor::view(Variant, rusanovMinusData);
          } else {
            neighKrnl.QgodNeighbor = runtime::init::QgodNeighbor::view(Variant, qGodNeighborData);
            neighKrnl.QcorrNeighbor = runtime::init::QcorrNeighbor::view(Variant, rusanovMinusNull);
          }
          neighKrnl.T = runtime::init::T::view(Variant, matTData);
          neighKrnl.Tinv = runtime::init::Tinv::view(Variant, matTinvData);
          neighKrnl.star(0) = runtime::init::star::view(Variant, 0, matATtildeData);
          if (boundaryProperties(cellInformation[cell].faceTypes[side]).usesFaceAlignedGhostState) {
            // the identity, in the layout it has as a tensor of its own
            neighKrnl.Tinv =
                runtime::init::identityT::view(Variant, init::identityT<Config>::Values);
          }
          neighKrnl.execute(Variant);

          if (cellInformation[cell].faceTypes[side] == FaceType::Dirichlet) {
            // the Dirichlet map is constant over the face, so it becomes part of
            // the local flux solver; what is left of the boundary condition is
            // the constant offset
            runtime::kernel::foldDirichlet foldKrnl;
            foldKrnl.AplusT =
                runtime::init::AplusT::view(Variant, localIntegration[cell].nApNm1[side]);
            foldKrnl.AminusT =
                runtime::init::AminusT::view(Variant, neighboringIntegration[cell].nAmNm1[side]);
            foldKrnl.Tinv = runtime::init::Tinv::view(Variant, matTinvData);
            foldKrnl.dirichletMap = runtime::init::dirichletMap::view(
                Variant, boundaryMapping[cell][side].dirichletMap);
            foldKrnl.execute(Variant);
          }

          if (cellInformation[cell].faceTypes[side] == FaceType::FreeSurfaceGravity) {
            // the free-surface-gravity map is constant over the face, so it becomes
            // part of the local flux solver; what is left of the boundary condition
            // is the displacement-driven offset
            runtime::kernel::foldFreeSurfaceGravity foldKrnl;
            foldKrnl.AplusT =
                runtime::init::AplusT::view(Variant, localIntegration[cell].nApNm1[side]);
            foldKrnl.AminusT =
                runtime::init::AminusT::view(Variant, neighboringIntegration[cell].nAmNm1[side]);
            foldKrnl.Tinv = runtime::init::Tinv::view(Variant, matTinvData);
            foldKrnl.execute(Variant);
          }
        }

        seissol::model::initializeSpecificLocalData<Config>(
            materialLocal, timeStepWidth, &localIntegration[cell].specific);

        seissol::model::initializeSpecificNeighborData<Config>(
            materialLocal, &neighboringIntegration[cell].specific);
      }
    }
  }
}

} // namespace seissol::initializer
