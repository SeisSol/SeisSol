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
#include "Equations/Impedance.h"      // IWYU pragma: keep
#include "Equations/ImpedanceBase.h"
#include "Equations/Setup.h" // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/Model/DynamicRuptureImpedance.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Model/CommonDatastructures.h"

#include <Eigen/Core>
#include <array>
#include <cassert>
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

/**
 * Copies an eigen3 matrix to a 2D yateto tensor; sparsely.
 */
template <typename T, typename S, int Dim1, int Dim2, size_t Dim1t>
void copyEigenToYateto(const Eigen::Matrix<T, Dim1, Dim2>& matrix,
                       yateto::CSCMatrixView<S, unsigned>& tensorView,
                       const std::array<size_t, Dim1t>& rowIdx) {
  // NOTE: shape(0) is the number of *logical* rows of the sparse matrix
  // (NumQuantities), not the number of stored ones. Only the row indices we
  // actually write have to be in range; they are ascending, so the last one
  // bounds them all.
  // (the dim parameters need to be int due to Eigen)
  static_assert(Dim1t == static_cast<size_t>(Dim1),
                "One target row index per row of the source matrix is required.");
  assert(rowIdx[Dim1t - 1] < tensorView.shape(0));
  // the target may be narrower than the source: the traction averaging matrices only carry the
  // components the frictional work is computed with, while the source also maps to the fluid
  // pressure. Writing past the pattern corrupts the neighboring entry.
  assert(tensorView.shape(1) <= static_cast<unsigned>(Dim2));

  tensorView.setZero();
  for (size_t row = 0; row < Dim1t; ++row) {
    for (size_t col = 0; col < tensorView.shape(1); ++col) {
      tensorView(rowIdx[row], col) = static_cast<S>(matrix(row, col));
    }
  }
}

/// The material of one side of a fault face, at the quadrature points of that
/// face.
///
/// A cell carries its material as samples; the fault wants it where its own
/// points are, and which parametrisation of the shared face applies depends on
/// the side: the face relation is 0 on the plus side and 1 on the minus side,
/// whose face orientation the canonical vertex order pins to zero. The
/// generated kernel folds the evaluation and the projection into one matrix per
/// side and face relation, so this runs
/// one small product per parameter the material declares. The samples are the
/// initialized ones the cell keeps, so nothing is fitted or derived a second
/// time here -- for a viscoelastic material that means the unrelaxed moduli
/// the fit produced, which is what the fault reads on the constant path too.
template <typename MaterialT>
void materialAtFaultPoints(const std::array<MaterialT, LTS::MaterialNodes>& samples,
                           std::uint8_t side,
                           std::uint8_t faceRelation,
                           const GlobalData& global,
                           std::array<MaterialT, seissol::dr::ImpedancePoints>& atPoints) {
  static_assert(tensor::materialAtFault::size() >= seissol::dr::ImpedancePoints,
                "The fault reads more points than the material is evaluated at.");

  alignas(Alignment) std::array<real, tensor::materialSamples::size()> sampled{};
  alignas(Alignment) std::array<real, tensor::materialAtFault::size()> atFault{};

  dynamicRupture::kernel::projectMaterialToFault krnl{};
  krnl.bindGlobals(global);
  krnl.materialSamples = sampled.data();
  krnl.materialAtFault = atFault.data();

  for (const auto& [name, member] : MaterialT::ParameterMap) {
    for (std::size_t node = 0; node < LTS::MaterialNodes; ++node) {
      // the material is one field, so every fused simulation sees the same
      // sample at a point
      for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
        sampled[sim + multisim::NumSimulations * node] = static_cast<real>(samples[node].*member);
      }
    }
    krnl.execute(side, faceRelation);
    for (std::size_t point = 0; point < atPoints.size(); ++point) {
      atPoints[point].*member = atFault[point];
    }
  }
}

/**
 * The "general" material case: impedance, eta and traction averaging matrices of a face whose
 * admittance is a full matrix, from the admittances of both sides.
 *
 * A template, so that the `if constexpr` below depends on MaterialT: every build instantiates it
 * with its own material, but the body is only compiled for the materials that take this path.
 * Their code generator gives the traction averaging matrices the full pattern, and their Riemann
 * problem couples the traction components. Everything else, isotropic elastic and viscoelastic
 * included, uses the scalar impedances.
 */
/// The scalar impedances of one point of a fault face, from the material on
/// either side of it there.
template <typename MaterialT>
void setIsotropicImpedance(seissol::dr::ImpedancesAndEta& impAndEta,
                           std::size_t point,
                           const MaterialT& plusMaterial,
                           const MaterialT& minusMaterial) {
  const double zp = plusMaterial.getDensity() * plusMaterial.getPWaveSpeed();
  const double zpNeig = minusMaterial.getDensity() * minusMaterial.getPWaveSpeed();
  const double zs = plusMaterial.getDensity() * plusMaterial.getSWaveSpeed();
  const double zsNeig = minusMaterial.getDensity() * minusMaterial.getSWaveSpeed();

  impAndEta.zp.set(point, zp);
  impAndEta.zpNeig.set(point, zpNeig);
  impAndEta.zs.set(point, zs);
  impAndEta.zsNeig.set(point, zsNeig);

  impAndEta.invZp.set(point, 1.0 / zp);
  impAndEta.invZpNeig.set(point, 1.0 / zpNeig);
  impAndEta.invZs.set(point, 1.0 / zs);
  impAndEta.invZsNeig.set(point, 1.0 / zsNeig);

  impAndEta.etaP.set(point, 1.0 / (1.0 / zp + 1.0 / zpNeig));
  impAndEta.invEtaS.set(point, 1.0 / zs + 1.0 / zsNeig);
  impAndEta.etaS.set(point, 1.0 / (1.0 / zs + 1.0 / zsNeig));
}

template <typename MaterialT>
void initializeFaultImpedance(const Fault& fault,
                              std::size_t meshFace,
                              std::size_t point,
                              const MaterialT& plusMaterial,
                              const MaterialT& minusMaterial,
                              seissol::dr::ImpedanceMatrices& impedanceMatrices,
                              DRGodunovData& godunovData,
                              seissol::dr::ImpedancesAndEta& impAndEta) {
  if constexpr (MaterialT::Type == seissol::model::MaterialType::Anisotropic ||
                MaterialT::Type == seissol::model::MaterialType::Poroelastic) {
    using ImpedanceCompute = seissol::model::ImpedanceCompute<MaterialT>;
    constexpr std::size_t N = ImpedanceCompute::Dim;
    // Zplus, Zminus and eta all share this dimension in the code generator
    static_assert(N == tensor::Zminus::Shape[0],
                  "The impedance tensors of the code generator do not match the material.");

    // the normal/tangent vectors are already normalized
    std::array<double, 36> bond{};
    seissol::model::getBondMatrix(fault.normal, fault.tangent1, fault.tangent2, bond);

    const auto plusLocal = seissol::model::getRotatedMaterialCoefficients(bond, plusMaterial);
    const auto minusLocal = seissol::model::getRotatedMaterialCoefficients(bond, minusMaterial);

    // Zplus/Zminus hold the *admittance* Y (traction -> velocity); eta is
    // (Y+ + Y-)^-1. For anisotropic materials Y is obtained in closed form
    // from the Christoffel matrix, which is exact also when qS1 and qS2 are
    // degenerate; for poroelasticity it comes from the Biot mass and stiffness
    // blocks in the same closed form.
    const auto faultImpedance =
        seissol::initializer::model::computeFaultImpedance(plusLocal, minusLocal);

    // The finite and consistency checks are a handful of flops per face and run in every
    // build: a material that is not positive definite produces NaN admittances right here,
    // and without the check the run only fails much later and somewhere else. Only the
    // self-adjointness and definiteness part costs an eigensolve, so that one stays behind
    // NDEBUG.
#ifdef NDEBUG
    constexpr bool CheckSelfAdjoint = false;
#else
    constexpr bool CheckSelfAdjoint = true;
#endif
    if (const auto violation =
            seissol::initializer::model::checkFaultImpedance(faultImpedance, CheckSelfAdjoint);
        violation.has_value()) {
      logError() << "Invalid dynamic rupture impedance at fault face" << meshFace << ":"
                 << violation.value();
    }

    const auto& impedanceMatrix = faultImpedance.admittancePlus;
    const auto& impedanceNeigMatrix = faultImpedance.admittanceMinus;
    const auto& etaMatrix = faultImpedance.eta;
    // the kernel contracts Q["kq"] * tractionMatrix["qp"], i.e. it applies
    // the transpose -- and b = eta * Y is not symmetric for a bimaterial
    // anisotropic interface (a few percent for realistic contrasts).
    const Eigen::Matrix<double, N, N> bMatrix = faultImpedance.bPlus.transpose();
    const Eigen::Matrix<double, N, N> bNeigMatrix = faultImpedance.bMinus.transpose();

    auto impedanceView = init::Zplus::view::create(impedanceMatrices.impedance.at(point));
    auto impedanceNeigView = init::Zminus::view::create(impedanceMatrices.impedanceNeig.at(point));
    auto etaView = init::eta::view::create(impedanceMatrices.eta.at(point));
    auto tractionPlusMatrix =
        init::tractionPlusMatrix::view::create(godunovData.tractionPlusMatrix);
    auto tractionMinusMatrix =
        init::tractionMinusMatrix::view::create(godunovData.tractionMinusMatrix);

    copyEigenToYateto(impedanceMatrix, impedanceView);
    copyEigenToYateto(impedanceNeigMatrix, impedanceNeigView);
    copyEigenToYateto(etaMatrix, etaView);
    // the rows of the traction averaging matrices; they have to match the sparsity pattern the
    // code generator builds for tractionPlusMatrix
    constexpr auto TractionRows = ImpedanceCompute::TractionIndices;
    copyEigenToYateto(bMatrix, tractionPlusMatrix, TractionRows);
    copyEigenToYateto(bNeigMatrix, tractionMinusMatrix, TractionRows);

    // reconstruction of the stress components outside of the Riemann problem; only needed by
    // the fault receiver output, which evaluates them on the plus side
    for (std::size_t col = 0; col < N; ++col) {
      for (std::size_t row = 0; row < 3; ++row) {
        impedanceMatrices.lateralStress.at(point)[col * 3 + row] =
            static_cast<real>(faultImpedance.lateralStressPlus(row, col));
      }
    }

    if constexpr (MaterialT::Type == seissol::model::MaterialType::Poroelastic) {
      // The solid frame is isotropic, so the shear rows of the Biot admittance decouple from
      // the fault-normal/fluid block and carry the same entry twice. A scalar impedance is
      // therefore exact here, and taking it from the admittance is what keeps the paths that
      // read ImpedancesAndEta -- the friction update, the slip accumulation and the receiver
      // output -- on the same Z_s = sqrt(mu * rho1) as the Riemann solver. Note that rho1 is
      // the statically condensed density, not the density of the solid grains.
      const double invZs = faultImpedance.admittancePlus(1, 1);
      const double invZsNeig = faultImpedance.admittanceMinus(1, 1);
      const double etaS = faultImpedance.eta(1, 1);

      impAndEta.zs.set(point, 1.0 / invZs);
      impAndEta.zsNeig.set(point, 1.0 / invZsNeig);
      impAndEta.invZs.set(point, invZs);
      impAndEta.invZsNeig.set(point, invZsNeig);
      impAndEta.etaS.set(point, etaS);
      impAndEta.invEtaS.set(point, 1.0 / etaS);
    }
  }
}

} // namespace

void initializeDynamicRuptureMatrices(const seissol::geometry::MeshReader& meshReader,
                                      LTS::Storage& ltsStorage,
                                      const LTS::Backmap& backmap,
                                      DynamicRupture::Storage& drStorage,
                                      const GlobalData& global) {
  real matTData[tensor::T::size()]{};
  real matTinvData[tensor::Tinv::size()]{};
  real matAPlusData[tensor::star::size(0)]{};
  real matAMinusData[tensor::star::size(0)]{};

  const auto& fault = meshReader.getFault();

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
      // Face relation 1 addresses the minus side at a zero face orientation index, which the
      // canonical vertex numbering guarantees on every interior face. Both sides of an MPI
      // split fault face therefore agree on it without exchanging anything.
      faceInformation[ltsFace].faceRelation = 1;
      faceInformation[ltsFace].plusSideOnThisRank = fault[meshFace].element >= 0;

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

      const auto getDofs = [&](const StoragePosition& position) -> real* {
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
          if (timeDerivative1 == nullptr &&
              cellInformation.ltsSetup.hasBuffer(BufferType::Derivatives)) {
            timeDerivative1 = ltsStorage.lookup<LTS::Derivatives>(position);
            timeDerivative1Device = ltsStorage.lookup<LTS::DerivativesDevice>(position);

            timeDofs1 = getDofs(position);
          }
          if (timeDerivative2 == nullptr &&
              cellInformation.ltsSetup.neighborBuffer(derivativesSide) == BufferType::Derivatives) {
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

      // The material at the points of the fault. Where it does not vary inside a
      // cell every point sees the cell's own material, so the two arrays hold
      // one entry and the loops below collapse.
      std::array<seissol::model::MaterialT, seissol::dr::ImpedancePoints> plusAtPoints{};
      std::array<seissol::model::MaterialT, seissol::dr::ImpedancePoints> minusAtPoints{};
      if constexpr (NodalMaterial) {
        // The samples of a cell that is not on this rank are reached through the
        // one that is: its face neighbour on the shared side, which is where the
        // ghost layer keeps them.
        const auto samplesOf = [&](const std::optional<StoragePosition>& own,
                                   const std::optional<StoragePosition>& other,
                                   std::uint8_t otherSide)
            -> const std::array<seissol::model::MaterialT, LTS::MaterialNodes>& {
          if (own.has_value()) {
            return ltsStorage.lookup<LTS::NodalMaterialData>(own.value());
          }
          const auto& secondary = ltsStorage.lookup<LTS::SecondaryInformation>(other.value());
          return ltsStorage.lookup<LTS::NodalMaterialData>(secondary.faceNeighbors[otherSide]);
        };

        // the plus side reads the face in its own parametrisation, the minus
        // side in the one the two sides agree on -- the same pair of arguments
        // the interpolation of the degrees of freedom uses
        materialAtFaultPoints(samplesOf(plusLtsId, minusLtsId, faceInformation[ltsFace].minusSide),
                              faceInformation[ltsFace].plusSide,
                              0,
                              global,
                              plusAtPoints);
        materialAtFaultPoints(samplesOf(minusLtsId, plusLtsId, faceInformation[ltsFace].plusSide),
                              faceInformation[ltsFace].minusSide,
                              faceInformation[ltsFace].faceRelation,
                              global,
                              minusAtPoints);
      } else {
        plusAtPoints[0] = *plusMaterial;
        minusAtPoints[0] = *minusMaterial;
      }

      // calculate Impedances Z and eta
      for (std::size_t point = 0; point < seissol::dr::ImpedancePoints; ++point) {
        setIsotropicImpedance(impAndEta[ltsFace], point, plusAtPoints[point], minusAtPoints[point]);
      }

      seissol::model::getTransposedCoefficientMatrix(*plusMaterial, 0, matAPlus);
      seissol::model::getTransposedCoefficientMatrix(*minusMaterial, 0, matAMinus);

      switch (plusMaterial->getMaterialType()) {
      case seissol::model::MaterialType::Anisotropic:
        [[fallthrough]];
      case seissol::model::MaterialType::Poroelastic: {
        for (std::size_t point = 0; point < seissol::dr::ImpedancePoints; ++point) {
          initializeFaultImpedance(fault[meshFace],
                                   meshFace,
                                   point,
                                   plusAtPoints[point],
                                   minusAtPoints[point],
                                   impedanceMatrices[ltsFace],
                                   godunovData[ltsFace],
                                   impAndEta[ltsFace]);
        }
        break;
      }
      default: {

        // NOTE: could be made `if constexpr`. However, that breaks ICC with a segfault.
        // So we don't do that, yet (until we drop ICC support at least).

        if (!::seissol::model::MaterialT::SupportsDR) {
          logError() << "The Dynamic Rupture mechanism does not work with the given material yet. "
                        "(built with:"
                     << ::seissol::model::MaterialT::Text << ")";
        }

        // the "fast" case, for isotropic elastic/viscoelastic. Does not need the extra impedance
        // matrices.

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
        break;
      }
      }

      /// Transpose matTinv.
      // Through the view rather than through a kernel, because TinvT is stored
      // in whichever layout the projections read it from -- packed to its
      // sparsity pattern where the build can take a packed operand -- and a
      // packed destination is not something the generated copy can write.
      // forall visits the entries the view actually stores, so the same line
      // fills a dense and a packed TinvT, and the entries a packed one leaves
      // out are the ones the rotation has no value for anyway.
      auto tinvT = init::TinvT::view::create(godunovData[ltsFace].dataTinvT);
      tinvT.forall(
          [&matTinv](const auto* entry, auto& value) { value = matTinv(entry[1], entry[0]); });

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
