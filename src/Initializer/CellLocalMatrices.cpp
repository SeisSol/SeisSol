/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Setup of SeisSol's cell local matrices.
 **/

#include "CellLocalMatrices.h"

#include <Initializer/BasicTypedefs.h>
#include <cassert>

#include "Equations/Setup.h" // IWYU pragma: keep
#include "Geometry/MeshTools.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/ParameterDB.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include "Parameters/ModelParameters.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"
#include <utils/logger.h>
#ifdef ACL_DEVICE
#include <device.h>
#endif

void setStarMatrix(
    const real* matAT, const real* matBT, const real* matCT, const real grad[3], real* starMatrix) {
  for (unsigned idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
    starMatrix[idx] = grad[0] * matAT[idx];
  }

  for (unsigned idx = 0; idx < seissol::tensor::star::size(1); ++idx) {
    starMatrix[idx] += grad[1] * matBT[idx];
  }

  for (unsigned idx = 0; idx < seissol::tensor::star::size(2); ++idx) {
    starMatrix[idx] += grad[2] * matCT[idx];
  }
}

void seissol::initializer::initializeCellLocalMatrices(
    const seissol::geometry::MeshReader& meshReader,
    LTSTree* ltsTree,
    LTS* lts,
    Lut* ltsLut,
    const TimeStepping& timeStepping,
    const parameters::ModelParameters& modelParameters) {
  const std::vector<Element>& elements = meshReader.getElements();
  const std::vector<Vertex>& vertices = meshReader.getVertices();

  static_assert(seissol::tensor::AplusT::Shape[0] == seissol::tensor::AminusT::Shape[0],
                "Shape mismatch for flux matrices");
  static_assert(seissol::tensor::AplusT::Shape[1] == seissol::tensor::AminusT::Shape[1],
                "Shape mismatch for flux matrices");

  assert(LayerMask(Ghost) == lts->material.mask);
  assert(LayerMask(Ghost) == lts->localIntegration.mask);
  assert(LayerMask(Ghost) == lts->neighboringIntegration.mask);

  const auto* cellInformationAll = ltsTree->var(lts->cellInformation);
  for (auto& layer : ltsTree->leaves(Ghost)) {
    auto* material = layer.var(lts->material);
    auto* localIntegration = layer.var(lts->localIntegration);
    auto* neighboringIntegration = layer.var(lts->neighboringIntegration);
    auto* cellInformation = layer.var(lts->cellInformation);
    auto* secondaryInformation = layer.var(lts->secondaryInformation);

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
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

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
        unsigned clusterId = secondaryInformation[cell].clusterId;
        auto timeStepWidth = timeStepping.globalCflTimeStepWidths[clusterId];
        unsigned meshId = secondaryInformation[cell].meshId;

        real x[4];
        real y[4];
        real z[4];
        real gradXi[3];
        real gradEta[3];
        real gradZeta[3];

        // Iterate over all 4 vertices of the tetrahedron
        for (unsigned vertex = 0; vertex < 4; ++vertex) {
          const VrtxCoords& coords = vertices[elements[meshId].vertices[vertex]].coords;
          x[vertex] = coords[0];
          y[vertex] = coords[1];
          z[vertex] = coords[2];
        }

        seissol::transformations::tetrahedronGlobalToReferenceJacobian(
            x, y, z, gradXi, gradEta, gradZeta);

        seissol::model::getTransposedCoefficientMatrix(material[cell].local, 0, matAT);
        seissol::model::getTransposedCoefficientMatrix(material[cell].local, 1, matBT);
        seissol::model::getTransposedCoefficientMatrix(material[cell].local, 2, matCT);
        setStarMatrix(
            matATData, matBTData, matCTData, gradXi, localIntegration[cell].starMatrices[0]);
        setStarMatrix(
            matATData, matBTData, matCTData, gradEta, localIntegration[cell].starMatrices[1]);
        setStarMatrix(
            matATData, matBTData, matCTData, gradZeta, localIntegration[cell].starMatrices[2]);

        double volume = MeshTools::volume(elements[meshId], vertices);

        for (unsigned side = 0; side < 4; ++side) {
          VrtxCoords normal;
          VrtxCoords tangent1;
          VrtxCoords tangent2;
          MeshTools::normalAndTangents(
              elements[meshId], side, vertices, normal, tangent1, tangent2);
          double surface = MeshTools::surface(normal);
          MeshTools::normalize(normal, normal);
          MeshTools::normalize(tangent1, tangent1);
          MeshTools::normalize(tangent2, tangent2);

          real nLocalData[6 * 6];
          seissol::model::getBondMatrix(normal, tangent1, tangent2, nLocalData);
          seissol::model::getTransposedGodunovState(
              seissol::model::getRotatedMaterialCoefficients(nLocalData, material[cell].local),
              seissol::model::getRotatedMaterialCoefficients(nLocalData,
                                                             material[cell].neighbor[side]),
              cellInformation[cell].faceTypes[side],
              qGodLocal,
              qGodNeighbor);
          seissol::model::getTransposedCoefficientMatrix(
              seissol::model::getRotatedMaterialCoefficients(nLocalData, material[cell].local),
              0,
              matATtilde);

          // Calculate transposed T instead
          seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, matT, matTinv);

          // Scale with |S_side|/|J| and multiply with -1 as the flux matrices
          // must be subtracted.
          real fluxScale = -2.0 * surface / (6.0 * volume);

          const auto isSpecialBC =
              [&secondaryInformation, &cellInformation, &cellInformationAll, cell](int side) {
                const auto hasDRFace = [](const CellLocalInformation& ci) {
                  bool hasAtLeastOneDRFace = false;
                  for (size_t i = 0; i < 4; ++i) {
                    if (ci.faceTypes[i] == FaceType::DynamicRupture) {
                      hasAtLeastOneDRFace = true;
                    }
                  }
                  return hasAtLeastOneDRFace;
                };
                const bool thisCellHasAtLeastOneDRFace = hasDRFace(cellInformation[cell]);
                const auto neighborID = secondaryInformation[cell].faceNeighborIds[side];
                const bool neighborBehindSideHasAtLeastOneDRFace =
                    hasDRFace(cellInformationAll[neighborID]);
                const bool adjacentDRFaceExists =
                    thisCellHasAtLeastOneDRFace || neighborBehindSideHasAtLeastOneDRFace;
                return (cellInformation[cell].faceTypes[side] == FaceType::Regular) &&
                       adjacentDRFaceExists;
              };

          const auto wavespeedLocal = material[cell].local.getMaxWaveSpeed();
          const auto wavespeedNeighbor = material[cell].neighbor[side].getMaxWaveSpeed();
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
              FaceType::FreeSurface, FaceType::FreeSurfaceGravity, FaceType::Analytical};

          const auto enforceGodunov = std::any_of(
              GodunovBoundaryConditions.begin(),
              GodunovBoundaryConditions.end(),
              [&](auto condition) { return condition == cellInformation[cell].faceTypes[side]; });

          const auto flux = enforceGodunov ? parameters::NumericalFlux::Godunov : fluxDefault;

          kernel::computeFluxSolverLocal localKrnl;
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
            material[cell].local, timeStepWidth, &localIntegration[cell].specific);

        seissol::model::initializeSpecificNeighborData(material[cell].local,
                                                       &neighboringIntegration[cell].specific);
      }
#ifdef _OPENMP
    }
#endif
  }
}

void surfaceAreaAndVolume(const seissol::geometry::MeshReader& meshReader,
                          unsigned meshId,
                          unsigned side,
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

void seissol::initializer::initializeBoundaryMappings(
    const seissol::geometry::MeshReader& meshReader,
    const EasiBoundary* easiBoundary,
    LTSTree* ltsTree,
    LTS* lts,
    Lut* ltsLut) {
  const std::vector<Element>& elements = meshReader.getElements();
  const std::vector<Vertex>& vertices = meshReader.getVertices();

  for (auto& layer : ltsTree->leaves(Ghost)) {
    auto* cellInformation = layer.var(lts->cellInformation);
    auto* boundary = layer.var(lts->boundaryMapping);
    auto* secondaryInformation = layer.var(lts->secondaryInformation);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      const auto& element = elements[secondaryInformation[cell].meshId];
      const double* coords[4];
      for (unsigned v = 0; v < 4; ++v) {
        coords[v] = vertices[element.vertices[v]].coords;
      }
      for (unsigned side = 0; side < 4; ++side) {
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
        for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
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
        real* matTData = boundary[cell][side].TData;
        real* matTinvData = boundary[cell][side].TinvData;
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
          easiBoundary->query(nodes, easiBoundaryMap, easiBoundaryConstant);
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

/**
 * Copies an eigen3 matrix to a 2D yateto tensor
 */
template <typename T, int Dim1, int Dim2>
void copyEigenToYateto(const Eigen::Matrix<T, Dim1, Dim2>& matrix,
                       yateto::DenseTensorView<2, T>& tensorView) {
  assert(tensorView.shape(0) == dim1);
  assert(tensorView.shape(1) == dim2);

  tensorView.setZero();
  for (size_t row = 0; row < Dim1; ++row) {
    for (size_t col = 0; col < Dim2; ++col) {
      tensorView(row, col) = matrix(row, col);
    }
  }
}

constexpr int N = tensor::Zminus::Shape[0];
Eigen::Matrix<real, N, N>
    extractMatrix(eigenvalues::Eigenpair<std::complex<double>,
                                         seissol::model::MaterialT::NumQuantities> eigenpair) {
#ifdef USE_POROELASTIC
  constexpr std::array<int, 4> TractionIndices = {0, 3, 5, 9};
  constexpr std::array<int, 4> VelocityIndices = {6, 7, 8, 10};
  constexpr std::array<int, 4> ColumnIndices = {0, 1, 2, 3};
#else
  constexpr std::array<int, 3> TractionIndices = {0, 3, 5};
  constexpr std::array<int, 3> VelocityIndices = {6, 7, 8};
  constexpr std::array<int, 3> ColumnIndices = {0, 1, 2};
#endif
  auto matrix = eigenpair.getVectorsAsMatrix();
  Eigen::Matrix<double, N, N> matRT = matrix(TractionIndices, ColumnIndices).real();
  Eigen::Matrix<double, N, N> matRTInv = matRT.inverse();
  Eigen::Matrix<double, N, N> matRU = matrix(VelocityIndices, ColumnIndices).real();
  Eigen::Matrix<double, N, N> matM = matRU * matRTInv;
  return matM.cast<real>();
};

void seissol::initializer::initializeDynamicRuptureMatrices(
    const seissol::geometry::MeshReader& meshReader,
    LTSTree* ltsTree,
    LTS* lts,
    Lut* ltsLut,
    LTSTree* dynRupTree,
    DynamicRupture* dynRup,
    unsigned* ltsFaceToMeshFace,
    const GlobalData& global,
    double etaHack) {
  real matTData[tensor::T::size()];
  real matTinvData[tensor::Tinv::size()];
  real matAPlusData[tensor::star::size(0)];
  real matAMinusData[tensor::star::size(0)];

  const std::vector<Fault>& fault = meshReader.getFault();
  const std::vector<Element>& elements = meshReader.getElements();
  CellDRMapping(*drMapping)[4] = ltsTree->var(lts->drMapping);
  CellDRMapping(*drMappingDevice)[4] = ltsTree->var(lts->drMappingDevice);
  CellMaterialData* material = ltsTree->var(lts->material);
  real** derivatives = ltsTree->var(lts->derivatives);
  real*(*faceNeighbors)[4] = ltsTree->var(lts->faceNeighbors);
  real** derivativesDevice = ltsTree->var(lts->derivativesDevice);
  real*(*faceNeighborsDevice)[4] = ltsTree->var(lts->faceNeighborsDevice);
  CellLocalInformation* cellInformation = ltsTree->var(lts->cellInformation);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (auto& layer : dynRupTree->leaves(Ghost)) {
    real** timeDerivativePlus = layer.var(dynRup->timeDerivativePlus);
    real** timeDerivativeMinus = layer.var(dynRup->timeDerivativeMinus);
    real** timeDerivativePlusDevice = layer.var(dynRup->timeDerivativePlusDevice);
    real** timeDerivativeMinusDevice = layer.var(dynRup->timeDerivativeMinusDevice);
    DRGodunovData* godunovData = layer.var(dynRup->godunovData);
    real(*imposedStatePlus)[tensor::QInterpolated::size()] =
        layer.var(dynRup->imposedStatePlus, AllocationPlace::Host);
    real(*imposedStateMinus)[tensor::QInterpolated::size()] =
        layer.var(dynRup->imposedStateMinus, AllocationPlace::Host);
    real(*fluxSolverPlus)[tensor::fluxSolver::size()] =
        layer.var(dynRup->fluxSolverPlus, AllocationPlace::Host);
    real(*fluxSolverMinus)[tensor::fluxSolver::size()] =
        layer.var(dynRup->fluxSolverMinus, AllocationPlace::Host);
    real(*imposedStatePlusDevice)[tensor::QInterpolated::size()] =
        layer.var(dynRup->imposedStatePlus, AllocationPlace::Device);
    real(*imposedStateMinusDevice)[tensor::QInterpolated::size()] =
        layer.var(dynRup->imposedStateMinus, AllocationPlace::Device);
    real(*fluxSolverPlusDevice)[tensor::fluxSolver::size()] =
        layer.var(dynRup->fluxSolverPlus, AllocationPlace::Device);
    real(*fluxSolverMinusDevice)[tensor::fluxSolver::size()] =
        layer.var(dynRup->fluxSolverMinus, AllocationPlace::Device);
    DRFaceInformation* faceInformation = layer.var(dynRup->faceInformation);
    seissol::model::IsotropicWaveSpeeds* waveSpeedsPlus = layer.var(dynRup->waveSpeedsPlus);
    seissol::model::IsotropicWaveSpeeds* waveSpeedsMinus = layer.var(dynRup->waveSpeedsMinus);
    seissol::dr::ImpedancesAndEta* impAndEta = layer.var(dynRup->impAndEta);
    seissol::dr::ImpedanceMatrices* impedanceMatrices = layer.var(dynRup->impedanceMatrices);

#ifdef _OPENMP
#pragma omp parallel for private(matTData, matTinvData, matAPlusData, matAMinusData)               \
    schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
      assert(fault[meshFace].element >= 0 || fault[meshFace].neighborElement >= 0);

      /// Face information
      faceInformation[ltsFace].meshFace = meshFace;
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
      int derivativesSide = 0;
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
      for (unsigned duplicate = 0; duplicate < Lut::MaxDuplicates; ++duplicate) {
        unsigned ltsId = ltsLut->ltsId(lts->cellInformation.mask, derivativesMeshId, duplicate);
        if (timeDerivative1 == nullptr && (cellInformation[ltsId].ltsSetup >> 9) % 2 == 1) {
          timeDerivative1 =
              derivatives[ltsLut->ltsId(lts->derivatives.mask, derivativesMeshId, duplicate)];
          timeDerivative1Device =
              derivativesDevice[ltsLut->ltsId(lts->derivatives.mask, derivativesMeshId, duplicate)];
        }
        if (timeDerivative2 == nullptr &&
            (cellInformation[ltsId].ltsSetup >> derivativesSide) % 2 == 1) {
          timeDerivative2 =
              faceNeighbors[ltsLut->ltsId(lts->faceNeighbors.mask, derivativesMeshId, duplicate)]
                           [derivativesSide];
          timeDerivative2Device = faceNeighborsDevice[ltsLut->ltsId(
              lts->faceNeighbors.mask, derivativesMeshId, duplicate)][derivativesSide];
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
      for (unsigned duplicate = 0; duplicate < Lut::MaxDuplicates; ++duplicate) {
        unsigned plusLtsId =
            (fault[meshFace].element >= 0)
                ? ltsLut->ltsId(lts->drMapping.mask, fault[meshFace].element, duplicate)
                : std::numeric_limits<unsigned>::max();
        unsigned minusLtsId =
            (fault[meshFace].neighborElement >= 0)
                ? ltsLut->ltsId(lts->drMapping.mask, fault[meshFace].neighborElement, duplicate)
                : std::numeric_limits<unsigned>::max();

        assert(duplicate != 0 || plusLtsId != std::numeric_limits<unsigned>::max() ||
               minusLtsId != std::numeric_limits<unsigned>::max());

        if (plusLtsId != std::numeric_limits<unsigned>::max()) {
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
          {
            CellDRMapping& mapping = drMapping[plusLtsId][faceInformation[ltsFace].plusSide];
            mapping.side = faceInformation[ltsFace].plusSide;
            mapping.faceRelation = 0;
            mapping.godunov = &imposedStatePlus[ltsFace][0];
            mapping.fluxSolver = &fluxSolverPlus[ltsFace][0];
            CellDRMapping& mappingDevice =
                drMappingDevice[plusLtsId][faceInformation[ltsFace].plusSide];
            mappingDevice.side = faceInformation[ltsFace].plusSide;
            mappingDevice.faceRelation = 0;
            mappingDevice.godunov = &imposedStatePlusDevice[ltsFace][0];
            mappingDevice.fluxSolver = &fluxSolverPlusDevice[ltsFace][0];
          }
        }
        if (minusLtsId != std::numeric_limits<unsigned>::max()) {
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
          {
            CellDRMapping& mapping = drMapping[minusLtsId][faceInformation[ltsFace].minusSide];
            mapping.side = faceInformation[ltsFace].minusSide;
            mapping.faceRelation = faceInformation[ltsFace].faceRelation;
            mapping.godunov = &imposedStateMinus[ltsFace][0];
            mapping.fluxSolver = &fluxSolverMinus[ltsFace][0];
            CellDRMapping& mappingDevice =
                drMappingDevice[minusLtsId][faceInformation[ltsFace].minusSide];
            mappingDevice.side = faceInformation[ltsFace].minusSide;
            mappingDevice.faceRelation = faceInformation[ltsFace].faceRelation;
            mappingDevice.godunov = &imposedStateMinusDevice[ltsFace][0];
            mappingDevice.fluxSolver = &fluxSolverMinusDevice[ltsFace][0];
          }
        }
      }

      /// Transformation matrix
      auto matT = init::T::view::create(matTData);
      auto matTinv = init::Tinv::view::create(matTinvData);
      seissol::model::getFaceRotationMatrix(fault[meshFace].normal,
                                            fault[meshFace].tangent1,
                                            fault[meshFace].tangent2,
                                            matT,
                                            matTinv);

      /// Materials
      seissol::model::Material* plusMaterial = nullptr;
      seissol::model::Material* minusMaterial = nullptr;
      unsigned plusLtsId = (fault[meshFace].element >= 0)
                               ? ltsLut->ltsId(lts->material.mask, fault[meshFace].element)
                               : std::numeric_limits<unsigned>::max();
      unsigned minusLtsId = (fault[meshFace].neighborElement >= 0)
                                ? ltsLut->ltsId(lts->material.mask, fault[meshFace].neighborElement)
                                : std::numeric_limits<unsigned>::max();

      assert(plusLtsId != std::numeric_limits<unsigned>::max() ||
             minusLtsId != std::numeric_limits<unsigned>::max());

      if (plusLtsId != std::numeric_limits<unsigned>::max()) {
        plusMaterial = &material[plusLtsId].local;
        minusMaterial = &material[plusLtsId].neighbor[faceInformation[ltsFace].plusSide];
      } else {
        assert(minusLtsId != std::numeric_limits<unsigned>::max());
        plusMaterial = &material[minusLtsId].neighbor[faceInformation[ltsFace].minusSide];
        minusMaterial = &material[minusLtsId].local;
      }

      /// Wave speeds and Coefficient Matrices
      auto matAPlus = init::star::view<0>::create(matAPlusData);
      auto matAMinus = init::star::view<0>::create(matAMinusData);

      waveSpeedsPlus[ltsFace].density = plusMaterial->rho;
      waveSpeedsMinus[ltsFace].density = minusMaterial->rho;
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
          etaHack / (1.0 / impAndEta[ltsFace].zp + 1.0 / impAndEta[ltsFace].zpNeig);
      impAndEta[ltsFace].invEtaS = 1.0 / impAndEta[ltsFace].zs + 1.0 / impAndEta[ltsFace].zsNeig;
      impAndEta[ltsFace].etaS =
          1.0 / (1.0 / impAndEta[ltsFace].zs + 1.0 / impAndEta[ltsFace].zsNeig);

      switch (plusMaterial->getMaterialType()) {
      case seissol::model::MaterialType::Acoustic: {
        logError() << "Dynamic Rupture does not work with an acoustic material.";
        break;
      }
      case seissol::model::MaterialType::Poroelastic: {
        // TODO (SW) Extract this into a function
        seissol::model::getTransposedCoefficientMatrix(
            *dynamic_cast<seissol::model::PoroElasticMaterial*>(plusMaterial), 0, matAPlus);
        seissol::model::getTransposedCoefficientMatrix(
            *dynamic_cast<seissol::model::PoroElasticMaterial*>(minusMaterial), 0, matAMinus);

        auto plusEigenpair = seissol::model::getEigenDecomposition(
            *dynamic_cast<seissol::model::PoroElasticMaterial*>(plusMaterial));
        auto minusEigenpair = seissol::model::getEigenDecomposition(
            *dynamic_cast<seissol::model::PoroElasticMaterial*>(minusMaterial));

        // The impedance matrices are diagonal in the (visco)elastic case, so we only store
        // the values Zp, Zs. In the poroelastic case, the fluid pressure and normal component
        // of the traction depend on each other, so we need a more complicated matrix structure.
        Eigen::Matrix<real, N, N> impedanceMatrix = extractMatrix(plusEigenpair);
        Eigen::Matrix<real, N, N> impedanceNeigMatrix = extractMatrix(minusEigenpair);
        Eigen::Matrix<real, N, N> etaMatrix = (impedanceMatrix + impedanceNeigMatrix).inverse();

        auto impedanceView = init::Zplus::view::create(impedanceMatrices[ltsFace].impedance);
        auto impedanceNeigView =
            init::Zminus::view::create(impedanceMatrices[ltsFace].impedanceNeig);
        auto etaView = init::eta::view::create(impedanceMatrices[ltsFace].eta);

        copyEigenToYateto(impedanceMatrix, impedanceView);
        copyEigenToYateto(impedanceNeigMatrix, impedanceNeigView);
        copyEigenToYateto(etaMatrix, etaView);

        break;
      }
      case seissol::model::MaterialType::Anisotropic: {
        logError() << "The Dynamic Rupture mechanism does not work with anisotropy yet.";
        // TODO(SW): Make DR work with anisotropy
        break;
      }
      case seissol::model::MaterialType::Elastic: {
        seissol::model::getTransposedCoefficientMatrix(
            *dynamic_cast<seissol::model::ElasticMaterial*>(plusMaterial), 0, matAPlus);
        seissol::model::getTransposedCoefficientMatrix(
            *dynamic_cast<seissol::model::ElasticMaterial*>(minusMaterial), 0, matAMinus);
        break;
      }
      case seissol::model::MaterialType::Viscoelastic: {
        seissol::model::getTransposedCoefficientMatrix(
            *dynamic_cast<seissol::model::ViscoElasticMaterial*>(plusMaterial), 0, matAPlus);
        seissol::model::getTransposedCoefficientMatrix(
            *dynamic_cast<seissol::model::ViscoElasticMaterial*>(minusMaterial), 0, matAMinus);
        break;
      }
      default: {
        logError() << "The Dynamic Rupture mechanism does not work with the given material yet.";
        break;
      }
      }
      /// Traction matrices for "average" traction
      auto tractionPlusMatrix =
          init::tractionPlusMatrix::view::create(godunovData[ltsFace].tractionPlusMatrix);
      auto tractionMinusMatrix =
          init::tractionMinusMatrix::view::create(godunovData[ltsFace].tractionMinusMatrix);
      double cZpP = plusMaterial->rho * waveSpeedsPlus[ltsFace].pWaveVelocity;
      double cZsP = plusMaterial->rho * waveSpeedsPlus[ltsFace].sWaveVelocity;
      double cZpM = minusMaterial->rho * waveSpeedsMinus[ltsFace].pWaveVelocity;
      double cZsM = minusMaterial->rho * waveSpeedsMinus[ltsFace].sWaveVelocity;
      double etaP = cZpP * cZpM / (cZpP + cZpM);
      double etaS = cZsP * cZsM / (cZsP + cZsM);

      tractionPlusMatrix.setZero();
      tractionPlusMatrix(0, 0) = etaP / cZpP;
      tractionPlusMatrix(3, 1) = etaS / cZsP;
      tractionPlusMatrix(5, 2) = etaS / cZsP;

      tractionMinusMatrix.setZero();
      tractionMinusMatrix(0, 0) = etaP / cZpM;
      tractionMinusMatrix(3, 1) = etaS / cZsM;
      tractionMinusMatrix(5, 2) = etaS / cZsM;

      /// Transpose matTinv
      dynamicRupture::kernel::transposeTinv ttKrnl;
      ttKrnl.Tinv = matTinvData;
      ttKrnl.TinvT = godunovData[ltsFace].TinvT;
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

      dynamicRupture::kernel::rotateFluxMatrix krnl;
      krnl.T = matTData;

      real(*fluxSolverPlusHost)[tensor::fluxSolver::size()] = layer.var(dynRup->fluxSolverPlus);
      real(*fluxSolverMinusHost)[tensor::fluxSolver::size()] = layer.var(dynRup->fluxSolverMinus);

      krnl.fluxSolver = fluxSolverPlusHost[ltsFace];
      krnl.fluxScaleDR = -2.0 * plusSurfaceArea / (6.0 * plusVolume);
      krnl.star(0) = matAPlusData;
      krnl.execute();

      krnl.fluxSolver = fluxSolverMinusHost[ltsFace];
      krnl.fluxScaleDR = 2.0 * minusSurfaceArea / (6.0 * minusVolume);
      krnl.star(0) = matAMinusData;
      krnl.execute();
    }

    layerLtsFaceToMeshFace += layer.getNumberOfCells();
  }
}

void seissol::initializer::copyCellMatricesToDevice(LTSTree* ltsTree,
                                                    LTS* lts,
                                                    LTSTree* dynRupTree,
                                                    DynamicRupture* dynRup,
                                                    LTSTree* boundaryTree,
                                                    Boundary* boundary) {
#ifdef ACL_DEVICE

#endif // ACL_DEVICE
}
