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
#include "Geometry/CellTransform.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BoundaryHelper.h"
#include "Initializer/BoundarySetup.h"
#include "Initializer/Model/CellFlux.h"
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

// The three directional star matrices share one shape and layout (they are
// clones of each other in the generator), so star(0) sizes all of them. A
// factored build only keeps star(0) in the generated code at all, for the
// kernels that set up the flux solvers.
void setStarMatrix(const real* matAT,
                   const real* matBT,
                   const real* matCT,
                   const std::array<double, Cell::Dim>& grad,
                   real* starMatrix) {
  for (std::size_t idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
    starMatrix[idx] = grad[0] * matAT[idx];
  }

  for (std::size_t idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
    starMatrix[idx] += grad[1] * matBT[idx];
  }

  for (std::size_t idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
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
                        double MaterialT::* member) {
    for (std::size_t node = 0; node < LTS::MaterialNodes; ++node) {
      // the material is one field, so every fused simulation sees the same
      // sample at a point
      for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
        samples[sim + multisim::NumSimulations * node] = static_cast<real>(source[node].*member);
      }
    }
  };
  // the simulation index is the leading dimension of the face values, and
  // every simulation holds the same material there
  static_assert(tensor::materialAtFace::size() >= multisim::NumSimulations * FluxFaceNodes,
                "The face values hold fewer nodes than the flux reads.");
  const auto scatter = [&](std::array<MaterialT, FluxFaceNodes>& target,
                           double MaterialT::* member) {
    for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
      target[node].*member = atFace[node * multisim::NumSimulations];
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

  // What a material holds beyond the parameters it binds is derived rather
  // than sampled, and only the bound parameters are interpolated to the face.
  // The rest comes from a sample as it is: right for what is the same at every
  // point, like the relaxation frequencies, and not interpolated for what is
  // not, like the theta of a viscoelastic material -- which the flux does not
  // read.
  own.fill(ownSamples[0]);
  neighbor.fill(neighborSamples[0]);

  for (const auto& [name, member] : MaterialT::ParameterMap) {
    fill(ownSamples, member);
    ownKrnl.execute(side);
    scatter(own, member);

    fill(neighborSamples, member);
    neighborKrnl.execute(faceRelation, neighborSide);
    scatter(neighbor, member);
  }
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
    auto* boundaryMapping = layer.var<LTS::BoundaryMapping>();
    auto* nodalMaterial = NodalMaterial ? layer.var<LTS::NodalMaterialData>() : nullptr;

#pragma omp parallel
    {
      real matATData[tensor::star::size(0)]{};
      real matATtildeData[tensor::star::size(0)]{};
      real matBTData[tensor::star::size(0)]{};
      real matCTData[tensor::star::size(0)]{};
      auto matAT = init::star::view<0>::create(matATData);
      // matAT with elastic parameters in local coordinate system, used for flux kernel
      auto matATtilde = init::star::view<0>::create(matATtildeData);
      auto matBT = init::star::view<0>::create(matBTData);
      auto matCT = init::star::view<0>::create(matCTData);

      real matTData[seissol::tensor::T::size()]{};
      real matTinvData[seissol::tensor::Tinv::size()]{};
      auto matT = init::T::view::create(matTData);
      auto matTinv = init::Tinv::view::create(matTinvData);

      // Where the ghost state is already rotated, the identity stands in for
      // Tinv. The flux solver reads it through Tinv's layout, which keeps only
      // the pattern of the rotation, so it is written in that layout rather
      // than handed over as the dense identityT.
      real identityTinvData[seissol::tensor::Tinv::size()]{};
      {
        const auto identity = init::identityT::view::create(init::identityT::Values);
        init::Tinv::view::create(identityTinvData).forall([&](const auto* entry, real& value) {
          value = identity.isInRange(entry[0], entry[1]) ? identity(entry[0], entry[1]) : 0;
        });
      }

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

        if constexpr (FactoredStar) {
          const std::array<std::array<double, Cell::Dim>, Cell::Dim> gradients{
              gradXi, gradEta, gradZeta};
          for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
            for (std::size_t component = 0; component < Cell::Dim; ++component) {
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
              if constexpr (NodalSource) {
                // the source term varies with the material just as the flux does
                const auto source = seissol::model::getSourceCoefficients(sampled[point]);
                for (std::size_t i = 0; i < source.size(); ++i) {
                  localIntegration[cell].sourceCoefficients[i][point] = source[i];
                }
                if constexpr (NodalSourceDeviation) {
                  // what this point asks for beyond the term the cell already
                  // carries, which is the part a solve done once per cell misses
                  const auto mean = seissol::model::getSourceCoefficients(materialLocal);
                  for (std::size_t i = 0; i < source.size(); ++i) {
                    localIntegration[cell].sourceDeviation[i][point] = source[i] - mean[i];
                  }
                }
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
          // Only the diagonal the Godunov state stores: a solver that folds the
          // relaxation into its quantities keeps the elastic rows alone, and
          // the views are that narrow.
          for (size_t i = 0; i < std::min(tensor::QgodLocal::Shape[0], tensor::QgodLocal::Shape[1]);
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

          if constexpr (NodalFlux) {
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

            for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
              std::array<double, FluxCoefficientCount> plus{};
              std::array<double, FluxCoefficientCount> minus{};
              fluxScalarsOfNode(ownAtFace[node],
                                neighborAtFace[node],
                                cellInformation[cell].faceTypes[side],
                                flux,
                                fluxScale,
                                plus,
                                minus);
              for (std::size_t c = 0; c < FluxCoefficientCount; ++c) {
                localIntegration[cell].fluxCoefficients[side][c][node] = plus[c];
                neighboringIntegration[cell].fluxCoefficients[side][c][node] = minus[c];
              }
            }
          }

          // The state the local flux applies: the Riemann problem's, or the central flux of the
          // Rusanov form, in the form the corrector of this build takes (toCorrectorForm). A
          // fault face has no Riemann problem of its own here, since the fault supplies its flux.
          const bool dynamicRupture =
              cellInformation[cell].faceTypes[side] == FaceType::DynamicRupture;
          real localStateData[tensor::QgodLocal::size()]{};
          if (!dynamicRupture) {
            std::copy_n(flux == parameters::NumericalFlux::Rusanov ? centralFluxData
                                                                   : qGodLocalData,
                        tensor::QgodLocal::size(),
                        localStateData);
          }
          auto localState = init::QgodLocal::view::create(localStateData);
          toCorrectorForm(localState);

          kernel::computeFluxSolverLocal localKrnl;
          localKrnl.fluxScale = fluxScale;
          localKrnl.AplusT = localIntegration[cell].nApNm1[side];
          localKrnl.QgodLocal = localStateData;
          localKrnl.QcorrLocal = flux == parameters::NumericalFlux::Rusanov && !dynamicRupture
                                     ? rusanovPlusData
                                     : rusanovPlusNull;
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
          if (boundaryProperties(cellInformation[cell].faceTypes[side]).usesFaceAlignedGhostState) {
            neighKrnl.Tinv = identityTinvData;
          }
          neighKrnl.execute();

          if (cellInformation[cell].faceTypes[side] == FaceType::Dirichlet) {
            // the Dirichlet map is constant over the face, so it becomes part of
            // the local flux solver; what is left of the boundary condition is
            // the constant offset
            kernel::foldDirichlet foldKrnl;
            foldKrnl.AplusT = localIntegration[cell].nApNm1[side];
            foldKrnl.AminusT = neighboringIntegration[cell].nAmNm1[side];
            foldKrnl.Tinv = matTinvData;
            foldKrnl.dirichletMap = boundaryMapping[cell][side].dirichletMap;
            foldKrnl.execute();
          }

          if (cellInformation[cell].faceTypes[side] == FaceType::FreeSurfaceGravity) {
            // the free-surface-gravity map is constant over the face, so it becomes
            // part of the local flux solver; what is left of the boundary condition
            // is the displacement-driven offset
            kernel::foldFreeSurfaceGravity foldKrnl;
            // fsgMap is a constant; only the pool holds it
            foldKrnl.bindGlobals(global);
            foldKrnl.AplusT = localIntegration[cell].nApNm1[side];
            foldKrnl.AminusT = neighboringIntegration[cell].nAmNm1[side];
            foldKrnl.Tinv = matTinvData;
            foldKrnl.execute();
          }
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
