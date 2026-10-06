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
#include "Common/Real.h"
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
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

namespace {

template <typename Cfg>
void setStarMatrix(const Real<Cfg>* matAT,
                   const Real<Cfg>* matBT,
                   const Real<Cfg>* matCT,
                   const std::array<double, Cell::Dim>& grad,
                   Real<Cfg>* starMatrix) {
  for (std::size_t idx = 0; idx < seissol::tensor::star<Cfg>::size(0); ++idx) {
    starMatrix[idx] = grad[0] * matAT[idx];
  }

  for (std::size_t idx = 0; idx < seissol::tensor::star<Cfg>::size(1); ++idx) {
    starMatrix[idx] += grad[1] * matBT[idx];
  }

  for (std::size_t idx = 0; idx < seissol::tensor::star<Cfg>::size(2); ++idx) {
    starMatrix[idx] += grad[2] * matCT[idx];
  }
}

/// The star matrices and flux solvers of the cells of a layer of the configuration `Cfg`.
template <typename Cfg>
void initializeCellLocalMatricesOfLayer(LTS::Layer& layer,
                                        const seissol::geometry::MeshReader& meshReader,
                                        LTS::Storage& ltsStorage,
                                        const ClusterLayout& clusterLayout,
                                        const parameters::ModelParameters& modelParameters) {
  using real = Real<Cfg>;
  const std::vector<Element>& elements = meshReader.getElements();
  const std::vector<Vertex>& vertices = meshReader.getVertices();
  constexpr auto Variant = configIdOf<Cfg>();

  static_assert(seissol::tensor::AplusT<Cfg>::Shape[0] == seissol::tensor::AminusT<Cfg>::Shape[0],
                "Shape mismatch for flux matrices");
  static_assert(seissol::tensor::AplusT<Cfg>::Shape[1] == seissol::tensor::AminusT<Cfg>::Shape[1],
                "Shape mismatch for flux matrices");

  auto* material = layer.var<LTS::Material>();
  auto* materialData = layer.var<LTS::MaterialData>(Cfg());
  auto* localIntegration = layer.var<LTS::LocalIntegration>(Cfg());
  auto* neighboringIntegration = layer.var<LTS::NeighboringIntegration>(Cfg());
  auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
  auto* boundaryMapping = layer.var<LTS::BoundaryMapping>(Cfg());

#pragma omp parallel
  {
    real matATData[tensor::star<Cfg>::size(0)]{};
    real matATtildeData[tensor::star<Cfg>::size(0)]{};
    real matBTData[tensor::star<Cfg>::size(1)]{};
    real matCTData[tensor::star<Cfg>::size(2)]{};
    auto matAT = init::star<Cfg>::template view<0>::create(matATData);
    // matAT with elastic parameters in local coordinate system, used for flux kernel
    auto matATtilde = init::star<Cfg>::template view<0>::create(matATtildeData);
    auto matBT = init::star<Cfg>::template view<0>::create(matBTData);
    auto matCT = init::star<Cfg>::template view<0>::create(matCTData);

    real matTData[seissol::tensor::T<Cfg>::size()]{};
    real matTinvData[seissol::tensor::Tinv<Cfg>::size()]{};
    auto matT = init::T<Cfg>::view::create(matTData);
    auto matTinv = init::Tinv<Cfg>::view::create(matTinvData);

    real qGodLocalData[tensor::QgodLocal<Cfg>::size()]{};
    real qGodNeighborData[tensor::QgodNeighbor<Cfg>::size()]{};
    auto qGodLocal = init::QgodLocal<Cfg>::view::create(qGodLocalData);
    auto qGodNeighbor = init::QgodNeighbor<Cfg>::view::create(qGodNeighborData);

    real rusanovPlusNull[tensor::QcorrLocal<Cfg>::size()]{};
    real rusanovMinusNull[tensor::QcorrNeighbor<Cfg>::size()]{};

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

      seissol::model::getTransposedCoefficientMatrix<Cfg>(materialLocal, 0, matAT);
      seissol::model::getTransposedCoefficientMatrix<Cfg>(materialLocal, 1, matBT);
      seissol::model::getTransposedCoefficientMatrix<Cfg>(materialLocal, 2, matCT);

      setStarMatrix<Cfg>(
          matATData, matBTData, matCTData, gradXi, localIntegration[cell].starMatrices[0]);
      setStarMatrix<Cfg>(
          matATData, matBTData, matCTData, gradEta, localIntegration[cell].starMatrices[1]);
      setStarMatrix<Cfg>(
          matATData, matBTData, matCTData, gradZeta, localIntegration[cell].starMatrices[2]);

      const double volume = MeshTools::volume(elements[meshId], vertices);

      for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
        CoordinateT normal{};
        CoordinateT tangent1{};
        CoordinateT tangent2{};
        MeshTools::normalAndTangents(elements[meshId], side, vertices, normal, tangent1, tangent2);
        const double surface = MeshTools::surface(normal);
        MeshTools::normalize(normal, normal);
        MeshTools::normalize(tangent1, tangent1);
        MeshTools::normalize(tangent2, tangent2);

        // the neighbor as a material of this cell, for the Riemann problem at their face; it may
        // compute in another configuration (checkConfigBoundaries admits the pairs that can)
        const auto neighborConfig = isInternalFaceType(cellInformation[cell].faceTypes[side])
                                        ? cellInformation[cell].neighborConfigIds[side]
                                        : configIdOf<Cfg>();
        const auto materialNeighbor = dispatchConfig(neighborConfig, [&](auto neighborCfg) {
          using NeighborMaterialT = model::MaterialOf<decltype(neighborCfg)>;
          if constexpr (model::CanNeighbor<model::MaterialOf<Cfg>, NeighborMaterialT>) {
            return model::neighborAs<model::MaterialOf<Cfg>>(
                dynamic_cast<const NeighborMaterialT&>(*material[cell].neighbor[side]));
          } else {
            logError() << "The materials" << model::MaterialOf<Cfg>::Text << "and"
                       << NeighborMaterialT::Text << "cannot be face neighbors.";
            return model::MaterialOf<Cfg>{};
          }
        });

        // Defines a rotation matrix for computing material properties in face-local coordinates
        // for anisotropy. It has no effect for isotropic materials.
        std::array<double, 36> nLocalData{};
        seissol::model::getBondMatrix(normal, tangent1, tangent2, nLocalData);
        seissol::model::getTransposedGodunovState(
            seissol::model::getRotatedMaterialCoefficients(nLocalData, materialLocal),
            seissol::model::getRotatedMaterialCoefficients(nLocalData, materialNeighbor),
            cellInformation[cell].faceTypes[side],
            qGodLocal,
            qGodNeighbor);
        seissol::model::getTransposedCoefficientMatrix<Cfg>(
            seissol::model::getRotatedMaterialCoefficients(nLocalData, materialLocal),
            0,
            matATtilde);

        // Calculate transposed T and Tinv instead
        seissol::model::getFaceRotationMatrix<Cfg>(normal, tangent1, tangent2, matT, matTinv);

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

        real centralFluxData[tensor::QgodLocal<Cfg>::size()]{};
        real rusanovPlusData[tensor::QcorrLocal<Cfg>::size()]{};
        real rusanovMinusData[tensor::QcorrNeighbor<Cfg>::size()]{};
        auto centralFluxView = init::QgodLocal<Cfg>::view::create(centralFluxData);
        auto rusanovPlusView = init::QcorrLocal<Cfg>::view::create(rusanovPlusData);
        auto rusanovMinusView = init::QcorrNeighbor<Cfg>::view::create(rusanovMinusData);
        // the diagonal ends where the stored block does: with the memory variables among the
        // unknowns, these matrices keep only the rows of the quantities of the Riemann problem
        for (size_t i = 0;
             i < std::min(tensor::QgodLocal<Cfg>::Shape[0], tensor::QgodLocal<Cfg>::Shape[1]);
             i++) {
          if (centralFluxView.isInRange(i, i)) {
            centralFluxView(i, i) = 0.5;
          }
          if (rusanovPlusView.isInRange(i, i)) {
            rusanovPlusView(i, i) = wavespeed * 0.5;
          }
          if (rusanovMinusView.isInRange(i, i)) {
            rusanovMinusView(i, i) = -wavespeed * 0.5;
          }
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
          neighKrnl.Tinv = runtime::init::identityT::view(Variant, init::identityT<Cfg>::Values);
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
          foldKrnl.dirichletMap =
              runtime::init::dirichletMap::view(Variant, boundaryMapping[cell][side].dirichletMap);
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

      seissol::model::initializeSpecificLocalData<Cfg>(
          materialLocal, timeStepWidth, &localIntegration[cell].specific);

      seissol::model::initializeSpecificNeighborData<Cfg>(materialLocal,
                                                          &neighboringIntegration[cell].specific);
    }
  }
}

} // namespace

void initializeCellLocalMatrices(const seissol::geometry::MeshReader& meshReader,
                                 LTS::Storage& ltsStorage,
                                 const ClusterLayout& clusterLayout,
                                 const parameters::ModelParameters& modelParameters) {
  assert(LayerMask(Ghost) == ltsStorage.info<LTS::Material>().mask);
  assert(LayerMask(Ghost) == ltsStorage.info<LTS::LocalIntegration>().mask);
  assert(LayerMask(Ghost) == ltsStorage.info<LTS::NeighboringIntegration>().mask);

  for (auto& layer : ltsStorage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      initializeCellLocalMatricesOfLayer<decltype(cfg)>(
          layer, meshReader, ltsStorage, clusterLayout, modelParameters);
    });
  }
}

} // namespace seissol::initializer
