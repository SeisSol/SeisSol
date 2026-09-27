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
#include "Initializer/BoundaryHelper.h"
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
#include <cstdint>
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

/// The material of both sides of a face, at the nodes of that face.
///
/// The cell reads its own samples through the face evaluation; the neighbour
/// sees the face with the other parametrisation, so its evaluation carries the
/// renumbering into this cell's ordering. Both give the same physical point at
/// the same index, which is what the flux of a node is built from.
template <typename MaterialT>
void faceMaterials(const GlobalData& global,
                   const std::array<MaterialT, LTS::MaterialNodes>& ownSamples,
                   const std::array<MaterialT, LTS::MaterialNodes>& neighborSamples,
                   std::uint8_t side,
                   std::uint8_t neighborSide,
                   std::uint8_t faceRelation,
                   std::array<MaterialT, FluxFaceNodes>& own,
                   std::array<MaterialT, FluxFaceNodes>& neighbor) {
  alignas(Alignment) std::array<real, tensor::materialSamples::size()> samples{};
  alignas(Alignment) std::array<real, tensor::materialAtFace::size()> atFace{};

  const auto fill = [&](const std::array<MaterialT, LTS::MaterialNodes>& source,
                        double MaterialT::*member) {
    for (std::size_t node = 0; node < LTS::MaterialNodes; ++node) {
      // the material is one field, so every fused simulation sees the same
      // sample at a point
      for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
        samples[sim + multisim::NumSimulations * node] = static_cast<real>(source[node].*member);
      }
    }
  };
  const auto scatter = [&](std::array<MaterialT, FluxFaceNodes>& target,
                           double MaterialT::*member) {
    for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
      target[node].*member = atFace[node];
    }
  };

  kernel::projectMaterialToFace ownKrnl{};
  ownKrnl.bindGlobals(global);
  ownKrnl.materialSamples = samples.data();
  ownKrnl.materialAtFace = atFace.data();

  kernel::projectMaterialToNeighborFace neighborKrnl{};
  neighborKrnl.bindGlobals(global);
  neighborKrnl.materialSamples = samples.data();
  neighborKrnl.materialAtFace = atFace.data();

  for (const auto& [name, member] : MaterialT::ParameterMap) {
    fill(ownSamples, member);
    ownKrnl.execute(side);
    scatter(own, member);

    fill(neighborSamples, member);
    neighborKrnl.execute(faceRelation, neighborSide);
    scatter(neighbor, member);
  }
}

/// The ten scalars the flux operator of one node is built from, in the
/// coordinates of the face.
///
/// The matrix form folds the rotation into what a face stores; here it stays
/// in the kernel, because rotated the operator no longer has ten degrees of
/// freedom but fifty-eight. What is stored is the operator as the Riemann
/// problem states it, read at the positions the generated table names.
///
/// The shape mirrors computeFluxSolverLocal and its neighbour exactly, down to
/// both sides contracting the Godunov state with the coefficient matrix of the
/// *local* material, and the scale the matrix form carries in AplusT riding on
/// the scalars instead.
template <typename MaterialT>
void fluxScalarsOfNode(const MaterialT& local,
                       const MaterialT& neighbor,
                       FaceType faceType,
                       parameters::NumericalFlux flux,
                       double fluxScale,
                       std::array<double, FluxCoefficientCount>& plus,
                       std::array<double, FluxCoefficientCount>& minus) {
  constexpr std::size_t N = MaterialT::NumQuantities;
  constexpr std::size_t Diagonal =
      std::min(tensor::QgodLocal::Shape[0], tensor::QgodLocal::Shape[1]);

  alignas(Alignment) std::array<real, tensor::QgodLocal::size()> godLocalData{};
  alignas(Alignment) std::array<real, tensor::QgodNeighbor::size()> godNeighborData{};
  auto godLocal = init::QgodLocal::view::create(godLocalData.data());
  auto godNeighbor = init::QgodNeighbor::view::create(godNeighborData.data());

  alignas(Alignment) std::array<real, tensor::star::size(0)> starData{};
  auto star = init::star::view<0>::create(starData.data());
  seissol::model::getTransposedCoefficientMatrix(local, 0, star);

  // the Riemann problem, or the central flux the Rusanov form uses instead
  double correction = 0.0;
  if (flux == parameters::NumericalFlux::Rusanov) {
    godLocal.setZero();
    godNeighbor.setZero();
    for (std::size_t i = 0; i < Diagonal; ++i) {
      godLocal(i, i) = 0.5;
      godNeighbor(i, i) = 0.5;
    }
    correction = std::max(local.getMaxWaveSpeed(), neighbor.getMaxWaveSpeed()) * 0.5;
  } else {
    seissol::model::getTransposedGodunovState(local, neighbor, faceType, godLocal, godNeighbor);
  }

  const auto read =
      [&](auto& godunov, double correctionSign, std::array<double, FluxCoefficientCount>& target) {
        for (std::size_t c = 0; c < FluxCoefficientCount; ++c) {
          const auto& source = generated::FluxCoefficientSources[c];
          double value = 0.0;
          for (std::size_t k = 0; k < N; ++k) {
            const double g = godunov.isInRange(source.row, k) ? godunov(source.row, k) : 0.0;
            const double a = star.isInRange(k, source.column) ? star(k, source.column) : 0.0;
            value += g * a;
          }
          // Qcorr is the diagonal the Rusanov form adds, and zero otherwise
          if (source.row == source.column && source.row < Diagonal) {
            value += correctionSign * correction;
          }
          target[c] = fluxScale * value;
        }
      };
  read(godLocal, 1.0, plus);
  read(godNeighbor, -1.0, minus);
}

void initializeCellLocalMatrices(const seissol::geometry::MeshReader& meshReader,
                                 LTS::Storage& ltsStorage,
                                 const ClusterLayout& clusterLayout,
                                 const parameters::ModelParameters& modelParameters,
                                 const GlobalData& global) {
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
    auto* nodalMaterial = NodalMaterial ? layer.var<LTS::NodalMaterialData>() : nullptr;

#pragma omp parallel
    {
      real matATData[tensor::star::size(0)]{};
      real matATtildeData[tensor::star::size(0)]{};
      real matBTData[tensor::star::size(1)]{};
      real matCTData[tensor::star::size(2)]{};
      auto matAT = init::star::view<0>::create(matATData);
      // matAT with elastic parameters in local coordinate system, used for flux kernel
      auto matATtilde = init::star::view<0>::create(matATtildeData);
      auto matBT = init::star::view<0>::create(matBTData);
      auto matCT = init::star::view<0>::create(matCTData);

      real matTData[seissol::tensor::T::size()]{};
      real matTinvData[seissol::tensor::Tinv::size()]{};
      auto matT = init::T::view::create(matTData);
      auto matTinv = init::Tinv::view::create(matTinvData);

      real qGodLocalData[tensor::QgodLocal::size()]{};
      real qGodNeighborData[tensor::QgodNeighbor::size()]{};
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

        if constexpr (FactoredStar) {
          const double* const gradients[3] = {gradXi, gradEta, gradZeta};
          for (std::size_t dim = 0; dim < 3; ++dim) {
            for (std::size_t component = 0; component < 3; ++component) {
              localIntegration[cell].referenceGradients[dim][component] = gradients[dim][component];
            }
          }
          if constexpr (NodalMaterial) {
            // one operator per sample point, so one coefficient per point
            const auto& sampled = layer.var<LTS::NodalMaterialData>()[cell];
            for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
              const auto coefficients = seissol::model::getStarCoefficients(sampled[point]);
              for (std::size_t i = 0; i < coefficients.size(); ++i) {
                localIntegration[cell].materialCoefficients[i][point] = coefficients[i];
              }
            }
          } else {
            const auto coefficients = seissol::model::getStarCoefficients(materialLocal);
            for (std::size_t i = 0; i < coefficients.size(); ++i) {
              localIntegration[cell].materialCoefficients[i][0] = coefficients[i];
            }
          }
        } else {
          seissol::model::getTransposedCoefficientMatrix(materialLocal, 0, matAT);
          seissol::model::getTransposedCoefficientMatrix(materialLocal, 1, matBT);
          seissol::model::getTransposedCoefficientMatrix(materialLocal, 2, matCT);

          setStarMatrix(
              matATData, matBTData, matCTData, gradXi, localIntegration[cell].starMatrices[0]);
          setStarMatrix(
              matATData, matBTData, matCTData, gradEta, localIntegration[cell].starMatrices[1]);
          setStarMatrix(
              matATData, matBTData, matCTData, gradZeta, localIntegration[cell].starMatrices[2]);
        }

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

          if constexpr (NodalMaterial) {
            // the operator at the nodes of the face, from the material of both
            // cells there; the rotation is the face's own and the same for both
            // sides, so it is kept once
            auto rotation = init::T::view::create(localIntegration[cell].faceRotation[side]);
            const auto source = init::T::view::create(matTData);
            rotation.setZero();
            for (std::size_t row = 0; row < tensor::T::Shape[0]; ++row) {
              for (std::size_t column = 0; column < tensor::T::Shape[1]; ++column) {
                if (rotation.isInRange(row, column)) {
                  rotation(row, column) = source(row, column);
                }
              }
            }

            const auto& ownSamples = nodalMaterial[cell];
            std::array<model::MaterialT, FluxFaceNodes> ownAtFace{};
            std::array<model::MaterialT, FluxFaceNodes> neighborAtFace{};
            if (isInternalFaceType(cellInformation[cell].faceTypes[side])) {
              const auto neighborPosition = secondaryInformation[cell].faceNeighbors[side];
              faceMaterials(global,
                            ownSamples,
                            ltsStorage.lookup<LTS::NodalMaterialData>(neighborPosition),
                            static_cast<std::uint8_t>(side),
                            cellInformation[cell].faceRelations[side][0],
                            cellInformation[cell].faceRelations[side][1],
                            ownAtFace,
                            neighborAtFace);
            } else {
              // a boundary face has no neighbour; the matrix form takes the
              // cell's own material for both sides and so does this
              faceMaterials(global,
                            ownSamples,
                            ownSamples,
                            static_cast<std::uint8_t>(side),
                            static_cast<std::uint8_t>(side),
                            0,
                            ownAtFace,
                            neighborAtFace);
              neighborAtFace = ownAtFace;
            }

            const bool dynamicRupture =
                cellInformation[cell].faceTypes[side] == FaceType::DynamicRupture;
            for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
              std::array<double, FluxCoefficientCount> plus{};
              std::array<double, FluxCoefficientCount> minus{};
              fluxScalarsOfNode(ownAtFace[node],
                                neighborAtFace[node],
                                cellInformation[cell].faceTypes[side],
                                flux,
                                dynamicRupture ? 0.0 : fluxScale,
                                plus,
                                minus);
              for (std::size_t c = 0; c < FluxCoefficientCount; ++c) {
                localIntegration[cell].fluxCoefficients[side][c][node] = plus[c];
                neighboringIntegration[cell].fluxCoefficients[side][c][node] = minus[c];
              }
            }
          }

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
