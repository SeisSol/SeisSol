// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "CellLocalMatrices.h"

#include "Equations/Setup.h" // IWYU pragma: keep
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshTools.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/ParameterDB.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include "Parameters/ModelParameters.h"
#include <Common/Constants.h>
#include <DynamicRupture/Typedefs.h>
#include <Equations/Datastructures.h> // IWYU pragma: keep
#include <GeneratedCode/init.h>
#include <Geometry/MeshDefinition.h>
#include <Geometry/MeshReader.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Backmap.h>
#include <Model/CommonDatastructures.h>
#include <Numerical/Eigenvalues.h>
#include <Solver/MultipleSimulations.h>
#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <limits>
#include <optional>
#include <utils/logger.h>
#include <vector>

namespace {

template <typename Cfg>
void setStarMatrix(const Real<Cfg>* matAT,
                   const Real<Cfg>* matBT,
                   const Real<Cfg>* matCT,
                   const double grad[3],
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

void surfaceAreaAndVolume(const seissol::geometry::MeshReader& meshReader,
                          std::size_t meshId,
                          int side,
                          double* surfaceArea,
                          double* volume) {
  const std::vector<Vertex>& vertices = meshReader.getVertices();
  const std::vector<Element>& elements = meshReader.getElements();

  VrtxCoords normal;
  VrtxCoords tangent1;
  VrtxCoords tangent2;
  MeshTools::normalAndTangents(elements[meshId], side, vertices, normal, tangent1, tangent2);

  *volume = MeshTools::volume(elements[meshId], vertices);
  *surfaceArea = MeshTools::surface(normal);
}

/**
 * Copies an eigen3 matrix to a 2D yateto tensor
 */
template <typename T, typename S, int Dim1, int Dim2>
void copyEigenToYateto(const Eigen::Matrix<T, Dim1, Dim2>& matrix,
                       yateto::DenseTensorView<2, S>& tensorView) {
  assert(tensorView.shape(0) == Dim1);
  assert(tensorView.shape(1) == Dim2);

  tensorView.setZero();
  for (size_t row = 0; row < Dim1; ++row) {
    for (size_t col = 0; col < Dim2; ++col) {
      tensorView(row, col) = static_cast<S>(matrix(row, col));
    }
  }
}

template <typename Cfg>
constexpr int N = tensor::Zminus<Cfg>::Shape[0];
template <typename T, typename Cfg>
Eigen::Matrix<T, N<Cfg>, N<Cfg>> extractMatrix(
    eigenvalues::Eigenpair<std::complex<double>, model::MaterialTT<Cfg>::NumQuantities> eigenpair) {
  std::vector<int> tractionIndices;
  std::vector<int> velocityIndices;
  std::vector<int> columnIndices;

  if constexpr (model::MaterialTT<Cfg>::Type == model::MaterialType::Poroelastic) {
    tractionIndices = {0, 3, 5, 9};
    velocityIndices = {6, 7, 8, 10};
    columnIndices = {0, 1, 2, 3};
  } else {
    tractionIndices = {0, 3, 5};
    velocityIndices = {6, 7, 8};
    columnIndices = {0, 1, 2};
  }

  auto matrix = eigenpair.getVectorsAsMatrix();
  const Eigen::Matrix<double, N<Cfg>, N<Cfg>> matRT = matrix(tractionIndices, columnIndices).real();
  const Eigen::Matrix<double, N<Cfg>, N<Cfg>> matRTInv = matRT.inverse();
  const Eigen::Matrix<double, N<Cfg>, N<Cfg>> matRU = matrix(velocityIndices, columnIndices).real();
  const Eigen::Matrix<double, N<Cfg>, N<Cfg>> matM = matRU * matRTInv;
  return matM.template cast<T>();
};

} // namespace

namespace seissol::initializer {

void initializeCellLocalMatrices(const seissol::geometry::MeshReader& meshReader,
                                 LTS::Storage& ltsStorage,
                                 const ClusterLayout& clusterLayout,
                                 const parameters::ModelParameters& modelParameters) {
  const std::vector<Element>& elements = meshReader.getElements();
  const std::vector<Vertex>& vertices = meshReader.getVertices();

  assert(LayerMask(Ghost) == ltsStorage.info<LTS::Material>().mask);
  assert(LayerMask(Ghost) == ltsStorage.info<LTS::LocalIntegration>().mask);
  assert(LayerMask(Ghost) == ltsStorage.info<LTS::NeighboringIntegration>().mask);

  const auto* cellInformationAll = ltsStorage.var<LTS::CellInformation>();
  for (auto& layer : ltsStorage.leaves(Ghost)) {
    layer.wrap([&](auto cfg) {
      using Cfg = decltype(cfg);
      using MaterialT = model::MaterialTT<Cfg>;

      static_assert(seissol::tensor::AplusT<Cfg>::Shape[0] ==
                        seissol::tensor::AminusT<Cfg>::Shape[0],
                    "Shape mismatch for flux matrices");
      static_assert(seissol::tensor::AplusT<Cfg>::Shape[1] ==
                        seissol::tensor::AminusT<Cfg>::Shape[1],
                    "Shape mismatch for flux matrices");

      // NOLINTNEXTLINE
      using real = Real<Cfg>;

      auto* material = layer.var<LTS::Material>();
      auto* materialData = layer.var<LTS::MaterialData>(cfg);
      auto* localIntegration = layer.var<LTS::LocalIntegration>(cfg);
      auto* neighboringIntegration = layer.var<LTS::NeighboringIntegration>(cfg);
      auto* cellInformation = layer.var<LTS::CellInformation>();
      auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        real matATData[tensor::star<Cfg>::size(0)];
        real matATtildeData[tensor::star<Cfg>::size(0)];
        real matBTData[tensor::star<Cfg>::size(1)];
        real matCTData[tensor::star<Cfg>::size(2)];
        auto matAT = init::star<Cfg>::template view<0>::create(matATData);
        // matAT with elastic parameters in local coordinate system, used for flux kernel
        auto matATtilde = init::star<Cfg>::template view<0>::create(matATtildeData);
        auto matBT = init::star<Cfg>::template view<0>::create(matBTData);
        auto matCT = init::star<Cfg>::template view<0>::create(matCTData);

        real matTData[seissol::tensor::T<Cfg>::size()];
        real matTinvData[seissol::tensor::Tinv<Cfg>::size()];
        auto matT = init::T<Cfg>::view::create(matTData);
        auto matTinv = init::Tinv<Cfg>::view::create(matTinvData);

        real qGodLocalData[tensor::QgodLocal<Cfg>::size()];
        real qGodNeighborData[tensor::QgodNeighbor<Cfg>::size()];
        auto qGodLocal = init::QgodLocal<Cfg>::view::create(qGodLocalData);
        auto qGodNeighbor = init::QgodNeighbor<Cfg>::view::create(qGodNeighborData);

        real rusanovPlusNull[tensor::QcorrLocal<Cfg>::size()]{};
        real rusanovMinusNull[tensor::QcorrNeighbor<Cfg>::size()]{};

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
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
            VrtxCoords normal;
            VrtxCoords tangent1;
            VrtxCoords tangent2;
            MeshTools::normalAndTangents(
                elements[meshId], side, vertices, normal, tangent1, tangent2);
            const double surface = MeshTools::surface(normal);
            MeshTools::normalize(normal, normal);
            MeshTools::normalize(tangent1, tangent1);
            MeshTools::normalize(tangent2, tangent2);

            double nLocalData[6 * 6];
            seissol::model::getBondMatrix(normal, tangent1, tangent2, nLocalData);
            seissol::model::getTransposedGodunovState<Cfg>(
                seissol::model::getRotatedMaterialCoefficients<Cfg>(nLocalData, materialLocal),
                seissol::model::getRotatedMaterialCoefficients<Cfg>(
                    nLocalData, *dynamic_cast<MaterialT*>(material[cell].neighbor[side])),
                cellInformation[cell].faceTypes[side],
                qGodLocal,
                qGodNeighbor);
            seissol::model::getTransposedCoefficientMatrix<Cfg>(
                seissol::model::getRotatedMaterialCoefficients<Cfg>(nLocalData, materialLocal),
                0,
                matATtilde);

            // Calculate transposed T instead
            seissol::model::getFaceRotationMatrix<Cfg>(normal, tangent1, tangent2, matT, matTinv);

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

            real centralFluxData[tensor::QgodLocal<Cfg>::size()]{};
            real rusanovPlusData[tensor::QcorrLocal<Cfg>::size()]{};
            real rusanovMinusData[tensor::QcorrNeighbor<Cfg>::size()]{};
            auto centralFluxView = init::QgodLocal<Cfg>::view::create(centralFluxData);
            auto rusanovPlusView = init::QcorrLocal<Cfg>::view::create(rusanovPlusData);
            auto rusanovMinusView = init::QcorrNeighbor<Cfg>::view::create(rusanovMinusData);
            for (size_t i = 0;
                 i < std::min(tensor::QgodLocal<Cfg>::Shape[0], tensor::QgodLocal<Cfg>::Shape[1]);
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

            kernel::computeFluxSolverLocal<Cfg> localKrnl;
            localKrnl.fluxScale = fluxScale;
            localKrnl.AplusT = localIntegration[cell].nApNm1[side];
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

            kernel::computeFluxSolverNeighbor<Cfg> neighKrnl;
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
              neighKrnl.Tinv = init::identityT<Cfg>::Values;
            }
            neighKrnl.execute();
          }

          seissol::model::initializeSpecificLocalData<Cfg>(
              materialLocal, timeStepWidth, &localIntegration[cell].specific);

          seissol::model::initializeSpecificNeighborData<Cfg>(
              materialLocal, &neighboringIntegration[cell].specific);
        }
      }
    });
  }
}

void initializeBoundaryMappings(const seissol::geometry::MeshReader& meshReader,
                                const EasiBoundaryT& easiBoundary,
                                LTS::Storage& ltsStorage) {
  const std::vector<Element>& elements = meshReader.getElements();
  const std::vector<Vertex>& vertices = meshReader.getVertices();

  for (auto& layer : ltsStorage.leaves(Ghost)) {
    layer.wrap([&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>;
      using MaterialT = model::MaterialTT<Cfg>;
      auto* cellInformation = layer.var<LTS::CellInformation>();
      auto* boundary = layer.var<LTS::BoundaryMapping>(cfg);
      auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
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
          real nodesReferenceData[nodal::tensor::nodes2D<Cfg>::Size];
          std::copy_n(nodal::init::nodes2D<Cfg>::Values,
                      nodal::tensor::nodes2D<Cfg>::Size,
                      nodesReferenceData);
          auto nodesReference = nodal::init::nodes2D<Cfg>::view::create(nodesReferenceData);
          auto* nodes = boundary[cell][side].nodes;
          assert(nodes != nullptr);
          auto offset = 0;
          for (std::size_t i = 0;
               i < nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisFunctionDimension];
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
          auto matT = init::T<Cfg>::view::create(matTData);
          auto matTinv = init::Tinv<Cfg>::view::create(matTinvData);

          VrtxCoords normal;
          VrtxCoords tangent1;
          VrtxCoords tangent2;
          MeshTools::normalAndTangents(element, side, vertices, normal, tangent1, tangent2);
          MeshTools::normalize(normal, normal);
          MeshTools::normalize(tangent1, tangent1);
          MeshTools::normalize(tangent2, tangent2);
          seissol::model::getFaceRotationMatrix<Cfg>(normal, tangent1, tangent2, matT, matTinv);

          // Evaluate easi boundary condition matrices if needed
          real* easiBoundaryMap = boundary[cell][side].easiBoundaryMap;
          real* easiBoundaryConstant = boundary[cell][side].easiBoundaryConstant;
          assert(easiBoundaryMap != nullptr);
          assert(easiBoundaryConstant != nullptr);
          if (cellInformation[cell].faceTypes[side] == FaceType::Dirichlet) {
            std::get<std::optional<EasiBoundary<Cfg>>>(easiBoundary)
                ->query(nodes, easiBoundaryMap, easiBoundaryConstant);
          } else {
            // Boundary should not be evaluated
            std::fill_n(easiBoundaryMap,
                        seissol::tensor::easiBoundaryMap<Cfg>::size(),
                        std::numeric_limits<real>::signaling_NaN());
            std::fill_n(easiBoundaryConstant,
                        seissol::tensor::easiBoundaryConstant<Cfg>::size(),
                        std::numeric_limits<real>::signaling_NaN());
          }
        }
      }
    });
  }
}

void initializeDynamicRuptureMatrices(const seissol::geometry::MeshReader& meshReader,
                                      LTS::Storage& ltsStorage,
                                      const LTS::Backmap& backmap,
                                      DynamicRupture::Storage& drStorage) {

  const auto& fault = meshReader.getFault();
  const auto& elements = meshReader.getElements();

  for (auto& layer : drStorage.leaves(Ghost)) {
    layer.wrap([&](auto cfg) {
      using Cfg = decltype(cfg);
      using MaterialT = model::MaterialTT<Cfg>;

      // NOLINTNEXTLINE
      using real = Real<Cfg>;

      real matTData[tensor::T<Cfg>::size()];
      real matTinvData[tensor::Tinv<Cfg>::size()];
      real matAPlusData[tensor::star<Cfg>::size(0)];
      real matAMinusData[tensor::star<Cfg>::size(0)];

      auto* timeDerivativePlus = layer.var<DynamicRupture::TimeDerivativePlus>(cfg);
      auto* timeDerivativeMinus = layer.var<DynamicRupture::TimeDerivativeMinus>(cfg);
      auto* timeDerivativePlusDevice = layer.var<DynamicRupture::TimeDerivativePlusDevice>(cfg);
      auto* timeDerivativeMinusDevice = layer.var<DynamicRupture::TimeDerivativeMinusDevice>(cfg);
      auto* godunovData = layer.var<DynamicRupture::GodunovData>(cfg);
      auto* imposedStatePlus =
          layer.var<DynamicRupture::ImposedStatePlus>(cfg, AllocationPlace::Host);
      auto* imposedStateMinus =
          layer.var<DynamicRupture::ImposedStateMinus>(cfg, AllocationPlace::Host);
      auto* fluxSolverPlus = layer.var<DynamicRupture::FluxSolverPlus>(cfg, AllocationPlace::Host);
      auto* fluxSolverMinus =
          layer.var<DynamicRupture::FluxSolverMinus>(cfg, AllocationPlace::Host);
      auto* imposedStatePlusDevice =
          layer.var<DynamicRupture::ImposedStatePlus>(cfg, AllocationPlace::Device);
      auto* imposedStateMinusDevice =
          layer.var<DynamicRupture::ImposedStateMinus>(cfg, AllocationPlace::Device);
      auto* fluxSolverPlusDevice =
          layer.var<DynamicRupture::FluxSolverPlus>(cfg, AllocationPlace::Device);
      auto* fluxSolverMinusDevice =
          layer.var<DynamicRupture::FluxSolverMinus>(cfg, AllocationPlace::Device);
      auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
      auto* waveSpeedsPlus = layer.var<DynamicRupture::WaveSpeedsPlus>();
      auto* waveSpeedsMinus = layer.var<DynamicRupture::WaveSpeedsMinus>();
      auto* impAndEta = layer.var<DynamicRupture::ImpAndEta>(cfg);
      auto* impedanceMatrices = layer.var<DynamicRupture::ImpedanceMatrices>(cfg);

#ifdef _OPENMP
#pragma omp parallel for private(matTData, matTinvData, matAPlusData, matAMinusData)               \
    schedule(static)
#endif
      for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        const std::size_t meshFace = faceInformation[ltsFace].meshFace;
        assert(fault[meshFace].element >= 0 || fault[meshFace].neighborElement >= 0);

        /// Face information
        // already set: faceInformation[ltsFace].meshFace = meshFace;
        faceInformation[ltsFace].plusSide = fault[meshFace].side;
        faceInformation[ltsFace].minusSide = fault[meshFace].neighborSide;
        if (fault[meshFace].element >= 0) {
          faceInformation[ltsFace].faceRelation =
              elements[fault[meshFace].element].sideOrientations[fault[meshFace].side] + 1;
          faceInformation[ltsFace].plusSideOnThisRank = true;
        } else {
          /// \todo check if this is correct
          faceInformation[ltsFace].faceRelation =
              elements[fault[meshFace].neighborElement]
                  .sideOrientations[fault[meshFace].neighborSide] +
              1;
          faceInformation[ltsFace].plusSideOnThisRank = false;
        }

        /// Look for time derivative mapping in all duplicates
        int derivativesMeshId = 0;
        unsigned derivativesSide = 0;
        if (fault[meshFace].element >= 0) {
          derivativesMeshId = fault[meshFace].element;
          derivativesSide = faceInformation[ltsFace].plusSide;
        } else {
          derivativesMeshId = fault[meshFace].neighborElement;
          derivativesSide = faceInformation[ltsFace].minusSide;
        }
        real* timeDerivative1 = nullptr;
        real* timeDerivative2 = nullptr;
        real* timeDerivative1Device = nullptr;
        real* timeDerivative2Device = nullptr;
        for (std::size_t duplicate = 0; duplicate < LTS::Backmap::MaxDuplicates; ++duplicate) {
          const auto positionOpt = backmap.getDup(derivativesMeshId, duplicate);
          if (positionOpt.has_value()) {
            const auto position = positionOpt.value();
            const auto& cellInformation = ltsStorage.lookup<LTS::CellInformation>(position);
            if (timeDerivative1 == nullptr && cellInformation.ltsSetup.hasDerivatives()) {
              timeDerivative1 = ltsStorage.lookup<LTS::Derivatives>(cfg, position);
              timeDerivative1Device = ltsStorage.lookup<LTS::DerivativesDevice>(cfg, position);
            }
            if (timeDerivative2 == nullptr &&
                cellInformation.ltsSetup.neighborHasDerivatives(derivativesSide)) {
              timeDerivative2 = reinterpret_cast<real*>(
                  ltsStorage.lookup<LTS::FaceNeighbors>(position)[derivativesSide]);
              timeDerivative2Device = reinterpret_cast<real*>(
                  ltsStorage.lookup<LTS::FaceNeighborsDevice>(position)[derivativesSide]);
            }
          }
        }

        assert(timeDerivative1 != nullptr && timeDerivative2 != nullptr);

        if (fault[meshFace].element >= 0) {
          timeDerivativePlus[ltsFace] = timeDerivative1;
          timeDerivativeMinus[ltsFace] = timeDerivative2;
          timeDerivativePlusDevice[ltsFace] = timeDerivative1Device;
          timeDerivativeMinusDevice[ltsFace] = timeDerivative2Device;
        } else {
          timeDerivativePlus[ltsFace] = timeDerivative2;
          timeDerivativeMinus[ltsFace] = timeDerivative1;
          timeDerivativePlusDevice[ltsFace] = timeDerivative2Device;
          timeDerivativeMinusDevice[ltsFace] = timeDerivative1Device;
        }

        assert(timeDerivativePlus[ltsFace] != nullptr && timeDerivativeMinus[ltsFace] != nullptr);

        /// DR mapping for elements
        for (std::size_t duplicate = 0; duplicate < LTS::Backmap::MaxDuplicates; ++duplicate) {
          const auto plusLtsId = (fault[meshFace].element >= 0)
                                     ? backmap.getDup(fault[meshFace].element, duplicate)
                                     : std::optional<StoragePosition>();
          const auto minusLtsId = (fault[meshFace].neighborElement >= 0)
                                      ? backmap.getDup(fault[meshFace].neighborElement, duplicate)
                                      : std::optional<StoragePosition>();

          assert(duplicate != 0 || plusLtsId.has_value() || minusLtsId.has_value());

          if (plusLtsId.has_value()) {
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
            {
              auto& mapping = ltsStorage.lookup<LTS::DRMapping>(
                  cfg, plusLtsId.value())[faceInformation[ltsFace].plusSide];
              mapping.side = faceInformation[ltsFace].plusSide;
              mapping.faceRelation = 0;
              mapping.godunov = &imposedStatePlus[ltsFace][0];
              mapping.fluxSolver = &fluxSolverPlus[ltsFace][0];
              auto& mappingDevice = ltsStorage.lookup<LTS::DRMappingDevice>(
                  cfg, plusLtsId.value())[faceInformation[ltsFace].plusSide];
              mappingDevice.side = faceInformation[ltsFace].plusSide;
              mappingDevice.faceRelation = 0;
              mappingDevice.godunov = &imposedStatePlusDevice[ltsFace][0];
              mappingDevice.fluxSolver = &fluxSolverPlusDevice[ltsFace][0];
            }
          }
          if (minusLtsId.has_value()) {
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
            {
              auto& mapping = ltsStorage.lookup<LTS::DRMapping>(
                  cfg, minusLtsId.value())[faceInformation[ltsFace].minusSide];
              mapping.side = faceInformation[ltsFace].minusSide;
              mapping.faceRelation = faceInformation[ltsFace].faceRelation;
              mapping.godunov = &imposedStateMinus[ltsFace][0];
              mapping.fluxSolver = &fluxSolverMinus[ltsFace][0];
              auto& mappingDevice = ltsStorage.lookup<LTS::DRMappingDevice>(
                  cfg, minusLtsId.value())[faceInformation[ltsFace].minusSide];
              mappingDevice.side = faceInformation[ltsFace].minusSide;
              mappingDevice.faceRelation = faceInformation[ltsFace].faceRelation;
              mappingDevice.godunov = &imposedStateMinusDevice[ltsFace][0];
              mappingDevice.fluxSolver = &fluxSolverMinusDevice[ltsFace][0];
            }
          }
        }

        /// Transformation matrix
        auto matT = init::T<Cfg>::view::create(matTData);
        auto matTinv = init::Tinv<Cfg>::view::create(matTinvData);
        seissol::model::getFaceRotationMatrix<Cfg>(fault[meshFace].normal,
                                                   fault[meshFace].tangent1,
                                                   fault[meshFace].tangent2,
                                                   matT,
                                                   matTinv);

        /// Materials
        const MaterialT* plusMaterial = nullptr;
        const MaterialT* minusMaterial = nullptr;
        const auto plusLtsId = (fault[meshFace].element >= 0)
                                   ? backmap.getDup(fault[meshFace].element, 0)
                                   : std::optional<StoragePosition>();
        const auto minusLtsId = (fault[meshFace].neighborElement >= 0)
                                    ? backmap.getDup(fault[meshFace].neighborElement, 0)
                                    : std::optional<StoragePosition>();

        assert(plusLtsId.has_value() || minusLtsId.has_value());

        if (plusLtsId.has_value()) {
          const auto& cellMaterialData = ltsStorage.lookup<LTS::Material>(plusLtsId.value());
          plusMaterial = dynamic_cast<const MaterialT*>(cellMaterialData.local);
          minusMaterial = dynamic_cast<const MaterialT*>(
              cellMaterialData.neighbor[faceInformation[ltsFace].plusSide]);
        } else {
          assert(minusLtsId.has_value());
          const auto& cellMaterialData = ltsStorage.lookup<LTS::Material>(minusLtsId.value());
          plusMaterial = dynamic_cast<const MaterialT*>(
              cellMaterialData.neighbor[faceInformation[ltsFace].minusSide]);
          minusMaterial = dynamic_cast<const MaterialT*>(cellMaterialData.local);
        }

        /// Wave speeds and Coefficient Matrices
        auto matAPlus = init::star<Cfg>::template view<0>::create(matAPlusData);
        auto matAMinus = init::star<Cfg>::template view<0>::create(matAMinusData);

        waveSpeedsPlus[ltsFace].density = plusMaterial->getDensity();
        waveSpeedsMinus[ltsFace].density = minusMaterial->getDensity();
        waveSpeedsPlus[ltsFace].pWaveVelocity = plusMaterial->getPWaveSpeed();
        waveSpeedsPlus[ltsFace].sWaveVelocity = plusMaterial->getSWaveSpeed();
        waveSpeedsMinus[ltsFace].pWaveVelocity = minusMaterial->getPWaveSpeed();
        waveSpeedsMinus[ltsFace].sWaveVelocity = minusMaterial->getSWaveSpeed();

        // calculate Impedances Z and eta
        impAndEta[ltsFace].zp =
            (waveSpeedsPlus[ltsFace].density * waveSpeedsPlus[ltsFace].pWaveVelocity);
        impAndEta[ltsFace].zpNeig =
            (waveSpeedsMinus[ltsFace].density * waveSpeedsMinus[ltsFace].pWaveVelocity);
        impAndEta[ltsFace].zs =
            (waveSpeedsPlus[ltsFace].density * waveSpeedsPlus[ltsFace].sWaveVelocity);
        impAndEta[ltsFace].zsNeig =
            (waveSpeedsMinus[ltsFace].density * waveSpeedsMinus[ltsFace].sWaveVelocity);

        impAndEta[ltsFace].invZp = 1 / impAndEta[ltsFace].zp;
        impAndEta[ltsFace].invZpNeig = 1 / impAndEta[ltsFace].zpNeig;
        impAndEta[ltsFace].invZs = 1 / impAndEta[ltsFace].zs;
        impAndEta[ltsFace].invZsNeig = 1 / impAndEta[ltsFace].zsNeig;

        impAndEta[ltsFace].etaP =
            1.0 / (1.0 / impAndEta[ltsFace].zp + 1.0 / impAndEta[ltsFace].zpNeig);
        impAndEta[ltsFace].invEtaS = 1.0 / impAndEta[ltsFace].zs + 1.0 / impAndEta[ltsFace].zsNeig;
        impAndEta[ltsFace].etaS =
            1.0 / (1.0 / impAndEta[ltsFace].zs + 1.0 / impAndEta[ltsFace].zsNeig);

        switch (plusMaterial->getMaterialType()) {
        case seissol::model::MaterialType::Poroelastic: {
          auto plusEigenpair = seissol::model::getEigenDecomposition<Cfg>(*plusMaterial);
          auto minusEigenpair = seissol::model::getEigenDecomposition<Cfg>(*minusMaterial);

          // The impedance matrices are diagonal in the (visco)elastic case, so we only store
          // the values Zp, Zs. In the poroelastic case, the fluid pressure and normal component
          // of the traction depend on each other, so we need a more complicated matrix structure.
          const Eigen::Matrix<double, N<Cfg>, N<Cfg>> impedanceMatrix =
              extractMatrix<double, Cfg>(plusEigenpair);
          const Eigen::Matrix<double, N<Cfg>, N<Cfg>> impedanceNeigMatrix =
              extractMatrix<double, Cfg>(minusEigenpair);
          const Eigen::Matrix<double, N<Cfg>, N<Cfg>> etaMatrix =
              (impedanceMatrix + impedanceNeigMatrix).inverse();

          auto impedanceView = init::Zplus<Cfg>::view::create(impedanceMatrices[ltsFace].impedance);
          auto impedanceNeigView =
              init::Zminus<Cfg>::view::create(impedanceMatrices[ltsFace].impedanceNeig);
          auto etaView = init::eta<Cfg>::view::create(impedanceMatrices[ltsFace].eta);

          copyEigenToYateto(impedanceMatrix, impedanceView);
          copyEigenToYateto(impedanceNeigMatrix, impedanceNeigView);
          copyEigenToYateto(etaMatrix, etaView);

          break;
        }
        default: {
          if constexpr (!MaterialT::SupportsDR) {
            logError()
                << "The Dynamic Rupture mechanism does not work with the given material yet. "
                   "(built with:"
                << MaterialT::Text << ")";
          }
          break;
        }
        }
        seissol::model::getTransposedCoefficientMatrix<Cfg>(*plusMaterial, 0, matAPlus);
        seissol::model::getTransposedCoefficientMatrix<Cfg>(*minusMaterial, 0, matAMinus);

        /// Traction matrices for "average" traction
        auto tractionPlusMatrix =
            init::tractionPlusMatrix<Cfg>::view::create(godunovData[ltsFace].tractionPlusMatrix);
        auto tractionMinusMatrix =
            init::tractionMinusMatrix<Cfg>::view::create(godunovData[ltsFace].tractionMinusMatrix);
        const double cZpP = plusMaterial->getDensity() * waveSpeedsPlus[ltsFace].pWaveVelocity;
        const double cZsP = plusMaterial->getDensity() * waveSpeedsPlus[ltsFace].sWaveVelocity;
        const double cZpM = minusMaterial->getDensity() * waveSpeedsMinus[ltsFace].pWaveVelocity;
        const double cZsM = minusMaterial->getDensity() * waveSpeedsMinus[ltsFace].sWaveVelocity;
        const double etaP = cZpP * cZpM / (cZpP + cZpM);
        const double etaS = cZsP * cZsM / (cZsP + cZsM);

        tractionPlusMatrix.setZero();
        tractionPlusMatrix(0, 0) = etaP / cZpP;
        tractionPlusMatrix(3, 1) = etaS / cZsP;
        tractionPlusMatrix(5, 2) = etaS / cZsP;

        tractionMinusMatrix.setZero();
        tractionMinusMatrix(0, 0) = etaP / cZpM;
        tractionMinusMatrix(3, 1) = etaS / cZsM;
        tractionMinusMatrix(5, 2) = etaS / cZsM;

        /// Transpose matTinv
        dynamicRupture::kernel::transposeTinv<Cfg> ttKrnl;
        ttKrnl.Tinv = matTinvData;
        ttKrnl.TinvT = godunovData[ltsFace].dataTinvT;
        ttKrnl.execute();

        double plusSurfaceArea = 0;
        double plusVolume = 0;
        double minusSurfaceArea = 0;
        double minusVolume = 0;
        double surfaceArea = 0;
        if (fault[meshFace].element >= 0) {
          surfaceAreaAndVolume(meshReader,
                               fault[meshFace].element,
                               fault[meshFace].side,
                               &plusSurfaceArea,
                               &plusVolume);
          surfaceArea = plusSurfaceArea;
        } else {
          /// Blow up solution on purpose if used by mistake
          plusSurfaceArea = 1.e99;
          plusVolume = 1.0;
        }
        if (fault[meshFace].neighborElement >= 0) {
          surfaceAreaAndVolume(meshReader,
                               fault[meshFace].neighborElement,
                               fault[meshFace].neighborSide,
                               &minusSurfaceArea,
                               &minusVolume);
          surfaceArea = minusSurfaceArea;
        } else {
          /// Blow up solution on purpose if used by mistake
          minusSurfaceArea = 1.e99;
          minusVolume = 1.0;
        }
        godunovData[ltsFace].doubledSurfaceArea = 2.0 * surfaceArea;

        dynamicRupture::kernel::rotateFluxMatrix<Cfg> krnl;
        krnl.T = matTData;

        real(*fluxSolverPlusHost)[tensor::fluxSolver<Cfg>::size()] =
            layer.var<DynamicRupture::FluxSolverPlus>(cfg);
        real(*fluxSolverMinusHost)[tensor::fluxSolver<Cfg>::size()] =
            layer.var<DynamicRupture::FluxSolverMinus>(cfg);

        krnl.fluxSolver = fluxSolverPlusHost[ltsFace];
        krnl.fluxScaleDR = -2.0 * plusSurfaceArea / (6.0 * plusVolume);
        krnl.star(0) = matAPlusData;
        krnl.execute();

        krnl.fluxSolver = fluxSolverMinusHost[ltsFace];
        krnl.fluxScaleDR = 2.0 * minusSurfaceArea / (6.0 * minusVolume);
        krnl.star(0) = matAMinusData;
        krnl.execute();
      }
    });
  }
}

} // namespace seissol::initializer
