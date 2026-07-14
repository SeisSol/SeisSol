// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "DynamicRuptureMatrices.h"

#include "DynamicRupture/Typedefs.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Model/CommonDatastructures.h"
#include "Numerical/Eigenvalues.h"

#include <Eigen/Core>
#include <array>
#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

namespace {

void surfaceAreaAndVolume(const seissol::geometry::MeshReader& meshReader,
                          std::size_t meshId,
                          std::int8_t side,
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
  // (the dim praameters need to be int due to Eigen)

  tensorView.setZero();
  for (size_t row = 0; row < Dim1; ++row) {
    for (size_t col = 0; col < Dim2; ++col) {
      tensorView(row, col) = static_cast<S>(matrix(row, col));
    }
  }
}

constexpr size_t N = tensor::Zminus::Shape[0];
template <typename T>
Eigen::Matrix<T, N, N>
    extractMatrix(eigenvalues::Eigenpair<std::complex<double>,
                                         seissol::model::MaterialT::NumQuantities> eigenpair) {
  std::vector<size_t> tractionIndices;
  std::vector<size_t> velocityIndices;
  std::vector<size_t> columnIndices;

  if constexpr (model::MaterialT::Type == model::MaterialType::Poroelastic) {
    tractionIndices = {0, 3, 5, 9};
    velocityIndices = {6, 7, 8, 10};
    columnIndices = {0, 1, 2, 3};
  } else {
    tractionIndices = {0, 3, 5};
    velocityIndices = {6, 7, 8};
    columnIndices = {0, 1, 2};
  }

  auto matrix = eigenpair.getVectorsAsMatrix();
  const Eigen::Matrix<double, N, N> matRT = matrix(tractionIndices, columnIndices).real();
  const Eigen::Matrix<double, N, N> matRTInv = matRT.inverse();
  const Eigen::Matrix<double, N, N> matRU = matrix(velocityIndices, columnIndices).real();
  const Eigen::Matrix<double, N, N> matM = matRU * matRTInv;
  return matM.cast<T>();
};

} // namespace

void initializeDynamicRuptureMatrices(const seissol::geometry::MeshReader& meshReader,
                                      LTS::Storage& ltsStorage,
                                      const LTS::Backmap& backmap,
                                      DynamicRupture::Storage& drStorage) {
  real matTData[tensor::T::size()]{};
  real matTinvData[tensor::Tinv::size()]{};
  real matAPlusData[tensor::star::size(0)]{};
  real matAMinusData[tensor::star::size(0)]{};

  const auto& fault = meshReader.getFault();
  const auto& elements = meshReader.getElements();

  for (auto& layer : drStorage.leaves(Ghost)) {
    auto* timeDofsPlus = layer.var<DynamicRupture::TimeDofsPlus>();
    auto* timeDofsMinus = layer.var<DynamicRupture::TimeDofsMinus>();
    auto* timeDerivativePlus = layer.var<DynamicRupture::TimeDerivativePlus>();
    auto* timeDerivativeMinus = layer.var<DynamicRupture::TimeDerivativeMinus>();
    auto* timeDerivativePlusDevice = layer.var<DynamicRupture::TimeDerivativePlusDevice>();
    auto* timeDerivativeMinusDevice = layer.var<DynamicRupture::TimeDerivativeMinusDevice>();
    auto* godunovData = layer.var<DynamicRupture::GodunovData>();
    auto* imposedStatePlus = layer.var<DynamicRupture::ImposedStatePlus>(AllocationPlace::Host);
    auto* imposedStateMinus = layer.var<DynamicRupture::ImposedStateMinus>(AllocationPlace::Host);
    auto* fluxSolverPlus = layer.var<DynamicRupture::FluxSolverPlus>(AllocationPlace::Host);
    auto* fluxSolverMinus = layer.var<DynamicRupture::FluxSolverMinus>(AllocationPlace::Host);
    auto* imposedStatePlusDevice =
        layer.var<DynamicRupture::ImposedStatePlus>(AllocationPlace::Device);
    auto* imposedStateMinusDevice =
        layer.var<DynamicRupture::ImposedStateMinus>(AllocationPlace::Device);
    auto* fluxSolverPlusDevice = layer.var<DynamicRupture::FluxSolverPlus>(AllocationPlace::Device);
    auto* fluxSolverMinusDevice =
        layer.var<DynamicRupture::FluxSolverMinus>(AllocationPlace::Device);
    auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
    auto* waveSpeedsPlus = layer.var<DynamicRupture::WaveSpeedsPlus>();
    auto* waveSpeedsMinus = layer.var<DynamicRupture::WaveSpeedsMinus>();
    auto* impAndEta = layer.var<DynamicRupture::ImpAndEta>();
    auto* impedanceMatrices = layer.var<DynamicRupture::ImpedanceMatrices>();

#pragma omp parallel for private(matTData, matTinvData, matAPlusData, matAMinusData)               \
    schedule(static)
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
      // TODO: change datatype after #1420
      int derivativesMeshId = 0;
      std::uint8_t derivativesSide = 0;
      if (fault[meshFace].element >= 0) {
        derivativesMeshId = fault[meshFace].element;
        derivativesSide = faceInformation[ltsFace].plusSide;
      } else {
        derivativesMeshId = fault[meshFace].neighborElement;
        derivativesSide = faceInformation[ltsFace].minusSide;
      }
      real* timeDofs1 = nullptr;
      real* timeDofs2 = nullptr;
      real* timeDerivative1 = nullptr;
      real* timeDerivative2 = nullptr;
      real* timeDerivative1Device = nullptr;
      real* timeDerivative2Device = nullptr;

      const auto getDofs = [&](const auto& position) -> real* {
        const auto halo = ltsStorage.getColorMap().argument(position.color).halo;
        if (halo == HaloType::Ghost) {
          return ltsStorage.lookup<LTS::DofsHalo>(position);
        } else {
          return ltsStorage.lookup<LTS::Dofs>(position);
        }
      };

      for (std::size_t duplicate = 0; duplicate < LTS::Backmap::MaxDuplicates; ++duplicate) {
        const auto positionOpt = backmap.getDup(derivativesMeshId, duplicate);
        if (positionOpt.has_value()) {
          const auto position = positionOpt.value();
          const auto& cellInformation = ltsStorage.lookup<LTS::CellInformation>(position);
          if (timeDerivative1 == nullptr && cellInformation.ltsSetup.hasDerivatives()) {
            timeDerivative1 = ltsStorage.lookup<LTS::Derivatives>(position);
            timeDerivative1Device = ltsStorage.lookup<LTS::DerivativesDevice>(position);

            timeDofs1 = getDofs(position);
          }
          if (timeDerivative2 == nullptr &&
              cellInformation.ltsSetup.neighborHasDerivatives(derivativesSide)) {
            timeDerivative2 = ltsStorage.lookup<LTS::FaceNeighbors>(position)[derivativesSide];
            timeDerivative2Device =
                ltsStorage.lookup<LTS::FaceNeighborsDevice>(position)[derivativesSide];

            const auto& secondaryInformation =
                ltsStorage.lookup<LTS::SecondaryInformation>(position);
            timeDofs2 = getDofs(secondaryInformation.faceNeighbors[derivativesSide]);
          }
        }
      }

      assert(timeDerivative1 != nullptr && timeDerivative2 != nullptr);

      if (fault[meshFace].element >= 0) {
        timeDofsPlus[ltsFace] = timeDofs1;
        timeDofsMinus[ltsFace] = timeDofs2;
        timeDerivativePlus[ltsFace] = timeDerivative1;
        timeDerivativeMinus[ltsFace] = timeDerivative2;
        timeDerivativePlusDevice[ltsFace] = timeDerivative1Device;
        timeDerivativeMinusDevice[ltsFace] = timeDerivative2Device;
      } else {
        timeDofsPlus[ltsFace] = timeDofs2;
        timeDofsMinus[ltsFace] = timeDofs1;
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

#pragma omp critical
          {
            CellDRMapping& mapping = ltsStorage.lookup<LTS::DRMapping>(
                plusLtsId.value())[faceInformation[ltsFace].plusSide];
            mapping.side = faceInformation[ltsFace].plusSide;
            mapping.faceRelation = 0;
            mapping.godunov = &imposedStatePlus[ltsFace][0];
            mapping.fluxSolver = &fluxSolverPlus[ltsFace][0];
            CellDRMapping& mappingDevice = ltsStorage.lookup<LTS::DRMappingDevice>(
                plusLtsId.value())[faceInformation[ltsFace].plusSide];
            mappingDevice.side = faceInformation[ltsFace].plusSide;
            mappingDevice.faceRelation = 0;
            mappingDevice.godunov = &imposedStatePlusDevice[ltsFace][0];
            mappingDevice.fluxSolver = &fluxSolverPlusDevice[ltsFace][0];
          }
        }
        if (minusLtsId.has_value()) {

#pragma omp critical
          {
            CellDRMapping& mapping = ltsStorage.lookup<LTS::DRMapping>(
                minusLtsId.value())[faceInformation[ltsFace].minusSide];
            mapping.side = faceInformation[ltsFace].minusSide;
            mapping.faceRelation = faceInformation[ltsFace].faceRelation;
            mapping.godunov = &imposedStateMinus[ltsFace][0];
            mapping.fluxSolver = &fluxSolverMinus[ltsFace][0];
            CellDRMapping& mappingDevice = ltsStorage.lookup<LTS::DRMappingDevice>(
                minusLtsId.value())[faceInformation[ltsFace].minusSide];
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
      const seissol::model::MaterialT* plusMaterial = nullptr;
      const seissol::model::MaterialT* minusMaterial = nullptr;
      const auto plusLtsId = (fault[meshFace].element >= 0)
                                 ? backmap.getDup(fault[meshFace].element, 0)
                                 : std::optional<StoragePosition>();
      const auto minusLtsId = (fault[meshFace].neighborElement >= 0)
                                  ? backmap.getDup(fault[meshFace].neighborElement, 0)
                                  : std::optional<StoragePosition>();

      assert(plusLtsId.has_value() || minusLtsId.has_value());

      if (plusLtsId.has_value()) {
        const auto& cellMaterialData = ltsStorage.lookup<LTS::Material>(plusLtsId.value());
        plusMaterial = dynamic_cast<seissol::model::MaterialT*>(cellMaterialData.local);
        minusMaterial = dynamic_cast<seissol::model::MaterialT*>(
            cellMaterialData.neighbor[faceInformation[ltsFace].plusSide]);
      } else {
        assert(minusLtsId.has_value());
        const auto& cellMaterialData = ltsStorage.lookup<LTS::Material>(minusLtsId.value());
        plusMaterial = dynamic_cast<seissol::model::MaterialT*>(
            cellMaterialData.neighbor[faceInformation[ltsFace].minusSide]);
        minusMaterial = dynamic_cast<seissol::model::MaterialT*>(cellMaterialData.local);
      }

      if (plusMaterial == nullptr || minusMaterial == nullptr) {
        logError() << "Materials on both sides of a fault face do not match.";
      }

      /// Wave speeds and Coefficient Matrices
      auto matAPlus = init::star::view<0>::create(matAPlusData);
      auto matAMinus = init::star::view<0>::create(matAMinusData);

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
        auto plusEigenpair = seissol::model::getEigenDecomposition(*plusMaterial);
        auto minusEigenpair = seissol::model::getEigenDecomposition(*minusMaterial);

        // The impedance matrices are diagonal in the (visco)elastic case, so we only store
        // the values Zp, Zs. In the poroelastic case, the fluid pressure and normal component
        // of the traction depend on each other, so we need a more complicated matrix structure.
        const Eigen::Matrix<double, N, N> impedanceMatrix = extractMatrix<double>(plusEigenpair);
        const Eigen::Matrix<double, N, N> impedanceNeigMatrix =
            extractMatrix<double>(minusEigenpair);
        const Eigen::Matrix<double, N, N> etaMatrix =
            (impedanceMatrix + impedanceNeigMatrix).inverse();

        auto impedanceView = init::Zplus::view::create(impedanceMatrices[ltsFace].impedance);
        auto impedanceNeigView =
            init::Zminus::view::create(impedanceMatrices[ltsFace].impedanceNeig);
        auto etaView = init::eta::view::create(impedanceMatrices[ltsFace].eta);

        copyEigenToYateto(impedanceMatrix, impedanceView);
        copyEigenToYateto(impedanceNeigMatrix, impedanceNeigView);
        copyEigenToYateto(etaMatrix, etaView);

        break;
      }
      default: {

        // NOTE: could be made `if constexpr`. However, that breaks ICC with a segfault.
        // So we don't do that, yet (until we drop ICC support at least).

        if (!model::MaterialT::SupportsDR) {
          logError() << "The Dynamic Rupture mechanism does not work with the given material yet. "
                        "(built with:"
                     << model::MaterialT::Text << ")";
        }
        break;
      }
      }
      seissol::model::getTransposedCoefficientMatrix(*plusMaterial, 0, matAPlus);
      seissol::model::getTransposedCoefficientMatrix(*minusMaterial, 0, matAMinus);

      /// Traction matrices for "average" traction
      auto tractionPlusMatrix =
          init::tractionPlusMatrix::view::create(godunovData[ltsFace].tractionPlusMatrix);
      auto tractionMinusMatrix =
          init::tractionMinusMatrix::view::create(godunovData[ltsFace].tractionMinusMatrix);
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
      dynamicRupture::kernel::transposeTinv ttKrnl;
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

      dynamicRupture::kernel::rotateFluxMatrix krnl;
      krnl.T = matTData;

      krnl.fluxSolver = fluxSolverPlus[ltsFace];
      krnl.fluxScaleDR = -2.0 * plusSurfaceArea / (6.0 * plusVolume);
      krnl.star(0) = matAPlusData;
      krnl.execute();

      krnl.fluxSolver = fluxSolverMinus[ltsFace];
      krnl.fluxScaleDR = 2.0 * minusSurfaceArea / (6.0 * minusVolume);
      krnl.star(0) = matAMinusData;
      krnl.execute();
    }
  }
}

} // namespace seissol::initializer
