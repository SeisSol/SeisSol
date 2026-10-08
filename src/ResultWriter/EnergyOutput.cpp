// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EnergyOutput.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "Equations/Energy.h" // IWYU pragma: keep
#include "Equations/EnergyBase.h"
#include "Equations/anisotropic/Model/Impedance.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "IO/Writer/File/RunFiles.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/PreProcessorMacros.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Modules/Modules.h"
#include "Monitoring/Unit.h"
#include "Numerical/Quadrature.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <mpi.h>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include "Common/Real.h"
#endif

namespace seissol::writer {

namespace {

template <typename Cfg>
std::array<Real<Cfg>, Cfg::NumSimulations>
    computeStaticWork(const Real<Cfg>* degreesOfFreedomPlus,
                      const Real<Cfg>* degreesOfFreedomMinus,
                      const DRFaceInformation& faceInfo,
                      const DRGodunovData<Cfg>& godunovData,
                      const Real<Cfg> slip[seissol::tensor::slipInterpolated<Cfg>::size()],
                      const GlobalData<Cfg>& global) {
  using real = Real<Cfg>;

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg> krnl;
  krnl.bindGlobals(global);

  alignas(PagesizeStack) real qInterpolatedPlus[tensor::QInterpolatedPlus<Cfg>::size()];
  alignas(PagesizeStack) real qInterpolatedMinus[tensor::QInterpolatedMinus<Cfg>::size()];
  alignas(Alignment) real tractionInterpolated[tensor::tractionInterpolated<Cfg>::size()];
  alignas(Alignment) real qPlus[tensor::Q<Cfg>::size()];
  alignas(Alignment) real qMinus[tensor::Q<Cfg>::size()];

  // needed to counter potential mis-alignment
  std::memcpy(qPlus, degreesOfFreedomPlus, sizeof(qPlus));
  std::memcpy(qMinus, degreesOfFreedomMinus, sizeof(qMinus));

  krnl.QInterpolated = qInterpolatedPlus;
  krnl.Q = qPlus;
  krnl.TinvT = godunovData.dataTinvT;
  krnl._prefetch.QInterpolated = qInterpolatedPlus;
  krnl.execute(faceInfo.plusSide, 0);

  krnl.QInterpolated = qInterpolatedMinus;
  krnl.Q = qMinus;
  krnl.TinvT = godunovData.dataTinvT;
  krnl._prefetch.QInterpolated = qInterpolatedMinus;
  krnl.execute(faceInfo.minusSide, faceInfo.faceRelation);

  constexpr auto Variant = configIdOf<Cfg>();
  runtime::dynamicRupture::kernel::computeTractionInterpolated trKrnl;
  trKrnl.tractionPlusMatrix =
      runtime::init::tractionPlusMatrix::view(Variant, godunovData.tractionPlusMatrix);
  trKrnl.tractionMinusMatrix =
      runtime::init::tractionMinusMatrix::view(Variant, godunovData.tractionMinusMatrix);
  trKrnl.QInterpolatedPlus = runtime::init::QInterpolatedPlus::view(Variant, qInterpolatedPlus);
  trKrnl.QInterpolatedMinus = runtime::init::QInterpolatedMinus::view(Variant, qInterpolatedMinus);
  trKrnl.tractionInterpolated =
      runtime::init::tractionInterpolated::view(Variant, tractionInterpolated);
  trKrnl.execute(Variant);

  alignas(Alignment) real staticFrictionalWork[tensor::staticFrictionalWork<Cfg>::size()]{};

  runtime::dynamicRupture::kernel::accumulateStaticFrictionalWork feKrnl;
  feKrnl.slipInterpolated = runtime::init::slipInterpolated::view(Variant, slip);
  feKrnl.tractionInterpolated =
      runtime::init::tractionInterpolated::view(Variant, tractionInterpolated);
  feKrnl.staticFrictionalWork =
      runtime::init::staticFrictionalWork::view(Variant, staticFrictionalWork);
  feKrnl.minusSurfaceArea = -0.5 * godunovData.doubledSurfaceArea;
  feKrnl.execute(Variant);

  std::array<real, Cfg::NumSimulations> frictionalWorkReturn{};
  std::copy_n(staticFrictionalWork, Cfg::NumSimulations, frictionalWorkReturn.begin());
  return frictionalWorkReturn;
}

// Energies that do not come from the material's EnergyCompute specialization.
// keep these here until we find a better place for them
//
// Note: plastic moment, frictional work and the momentum triple are printed by
// bespoke blocks in printEnergies (they need an equivalent magnitude, a
// three-component grouping, or a ratio against a non-adjacent quantity), so
// their descriptors carry no label.

constexpr std::string_view PlasticMoment = "plastic_moment";
constexpr std::string_view GravitationalPotentialEnergy = "gravitational_potential_energy";
constexpr std::string_view SeismicMoment = "seismic_moment";
constexpr std::string_view TotalFrictionalWork = "total_frictional_work";
constexpr std::string_view StaticFrictionalWork = "static_frictional_work";
constexpr std::string_view Potency = "potency";

constexpr std::array GlobalEnergies{
    model::EnergyDescriptor{PlasticMoment, model::EnergyUnit::Moment, {}, {}, {}},
    model::EnergyDescriptor{GravitationalPotentialEnergy,
                            model::EnergyUnit::Energy,
                            "gravitational",
                            "Gravitational potential energy:",
                            {}},
    model::EnergyDescriptor{SeismicMoment, model::EnergyUnit::Moment, {}, {}, {}},
    model::EnergyDescriptor{TotalFrictionalWork, model::EnergyUnit::Energy, {}, {}, {}},
    model::EnergyDescriptor{StaticFrictionalWork, model::EnergyUnit::Energy, {}, {}, {}},
    model::EnergyDescriptor{Potency, model::EnergyUnit::Scalar, {}, {}, {}},
};
static_assert(model::detail::descriptorsWellFormed(GlobalEnergies),
              "energy descriptors must be named, unique, and grouped consistently");

constexpr std::array MomentumComponents{
    std::string_view{"momentumX"}, std::string_view{"momentumY"}, std::string_view{"momentumZ"}};

/// Adds the frictional work, the seismic moment and the potency of the faces of `layer`, which
/// compute in the configuration `Cfg`, and lowers the times since the slip rate fell below its
/// threshold to the ones of these faces.
template <typename Cfg>
void addFaultEnergies(const DynamicRupture::Layer& layer,
                      const GlobalData<Cfg>& global,
                      EnergiesStorage& energies,
                      std::vector<double>& minTimeSinceSlipRateBelowThreshold) {
  using real = Real<Cfg>;
  using MaterialT = model::MaterialOf<Cfg>;
  constexpr auto SimCount = Cfg::NumSimulations;

  double totalFrictionalWork[SimCount]{};
  double staticFrictionalWork[SimCount]{};
  double seismicMoment[SimCount]{};
  double potency[SimCount]{};

  real* const* timeDofsPlus = layer.var<DynamicRupture::TimeDerivativePlus>(Cfg());
  real* const* timeDofsMinus = layer.var<DynamicRupture::TimeDerivativeMinus>(Cfg());

  const auto* godunovData = layer.var<DynamicRupture::GodunovData>(Cfg());
  const auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
  const auto* drEnergyOutput = layer.var<DynamicRupture::DREnergyOutputVar>(Cfg());
  const auto* waveSpeedsPlus = layer.var<DynamicRupture::WaveSpeedsPlus>();
  const auto* waveSpeedsMinus = layer.var<DynamicRupture::WaveSpeedsMinus>();
  const auto* impedanceMatrices = layer.var<DynamicRupture::ImpedanceMatrices>(Cfg());
  const auto layerSize = layer.size();

#if !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(+ : totalFrictionalWork[ : SimCount],                           \
                                       staticFrictionalWork[ : SimCount],                          \
                                       seismicMoment[ : SimCount],                                 \
                                       potency[ : SimCount])
#endif
  for (std::size_t i = 0; i < layerSize; ++i) {
    if (faceInformation[i].plusSideOnThisRank) {
      const auto staticFrictionalWorkIncrease = computeStaticWork<Cfg>(timeDofsPlus[i],
                                                                       timeDofsMinus[i],
                                                                       faceInformation[i],
                                                                       godunovData[i],
                                                                       drEnergyOutput[i].slip,
                                                                       global);

#pragma omp simd
      for (size_t sim = 0; sim < SimCount; sim++) {
        staticFrictionalWork[sim] += staticFrictionalWorkIncrease[sim];
        for (std::size_t j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints<Cfg>; ++j) {
          totalFrictionalWork[sim] += drEnergyOutput[i].frictionalEnergy[j * SimCount + sim];
        }

        const double areaWeight = godunovData[i].doubledSurfaceArea;
        double potencyIncrease = 0.0;
        double momentIncrease = 0.0;

        if constexpr (MaterialT::Type == model::MaterialType::Anisotropic) {
          // The modulus turning potency into moment is d^T Gamma d with the fault-local
          // Christoffel matrix Gamma and the unit slip direction d, i.e. the contraction of the
          // moment tensor C_ijkl (n_k d_l + n_l d_k) / 2 with the source geometry. Both the
          // orientation of the fault and the rake enter, so the modulus varies from point to
          // point and cannot be pulled out of the quadrature sum.

          static_assert(MaterialT::Type != model::MaterialType::Anisotropic ||
                        (tensor::Zplus<Cfg>::size() == 9 && tensor::Zminus<Cfg>::size() == 9));

          const auto admittance = [](const real* data) {
            return Eigen::Map<const Eigen::Matrix<real, 3, 3>>(data).template cast<double>();
          };
          using AnisotropicImpedance = model::ImpedanceCompute<model::AnisotropicMaterial>;
          const auto gammaPlus = AnisotropicImpedance::christoffelFromAdmittance(
              admittance(impedanceMatrices[i].impedance), waveSpeedsPlus[i].density);
          const auto gammaMinus = AnisotropicImpedance::christoffelFromAdmittance(
              admittance(impedanceMatrices[i].impedanceNeig), waveSpeedsMinus[i].density);

          const auto* slip =
              reinterpret_cast<const real(*)[seissol::dr::misc::NumPaddedPoints<Cfg>]>(
                  drEnergyOutput[i].slip);

          for (std::size_t k = 0; k < seissol::dr::misc::NumBoundaryGaussPoints<Cfg>; ++k) {
            const auto index = k * SimCount + sim;

            // the rake is taken from the net slip; it is the instantaneous one only as long as
            // the slip direction does not turn during rupture
            const double slipStrike = slip[1][index];
            const double slipDip = slip[2][index];
            const double magnitude = std::sqrt(slipStrike * slipStrike + slipDip * slipDip);
            const double d1 = magnitude > 0 ? slipStrike / magnitude : 1.0;
            const double d2 = magnitude > 0 ? slipDip / magnitude : 0.0;

            const auto project = [d1, d2](const Eigen::Matrix3d& gamma) {
              return gamma(1, 1) * d1 * d1 + (gamma(1, 2) + gamma(2, 1)) * d1 * d2 +
                     gamma(2, 2) * d2 * d2;
            };
            const double muPlus = project(gammaPlus);
            const double muMinus = project(gammaMinus);

            const double slipIncrease =
                drEnergyOutput[i].accumulatedSlip[index] * init::quadweights<Cfg>::Values[k];
            potencyIncrease += slipIncrease;
            momentIncrease += slipIncrease * 2.0 * muPlus * muMinus / (muPlus + muMinus);
          }
          potencyIncrease *= areaWeight;
          momentIncrease *= areaWeight;
        } else {
          // rho * cs^2 is the shear modulus of the frame for every material with an isotropic
          // one, poroelasticity included -- there the fluid carries no shear
          const double muPlus = waveSpeedsPlus[i].density * waveSpeedsPlus[i].sWaveVelocity *
                                waveSpeedsPlus[i].sWaveVelocity;
          const double muMinus = waveSpeedsMinus[i].density * waveSpeedsMinus[i].sWaveVelocity *
                                 waveSpeedsMinus[i].sWaveVelocity;
          const double mu = 2.0 * muPlus * muMinus / (muPlus + muMinus);
          for (std::size_t k = 0; k < seissol::dr::misc::NumBoundaryGaussPoints<Cfg>; ++k) {
            potencyIncrease += drEnergyOutput[i].accumulatedSlip[k * SimCount + sim] *
                               init::quadweights<Cfg>::Values[k];
          }
          potencyIncrease *= areaWeight;
          momentIncrease = potencyIncrease * mu;
        }

        potency[sim] += potencyIncrease;
        seismicMoment[sim] += momentIncrease;
      }
    }
  }

  double localMin[SimCount]{};
  for (std::size_t sim = 0; sim < SimCount; ++sim) {
    localMin[sim] = std::numeric_limits<double>::max();
  }

#if !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(min : localMin[ : SimCount]) default(none)                      \
    shared(layerSize, drEnergyOutput, faceInformation, SimCount)
#endif
  for (std::size_t i = 0; i < layerSize; ++i) {
    if (faceInformation[i].plusSideOnThisRank) {

#pragma omp simd
      for (size_t sim = 0; sim < SimCount; sim++) {
        for (std::size_t j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints<Cfg>; ++j) {
          localMin[sim] = std::min(
              static_cast<double>(
                  drEnergyOutput[i]
                      .timeSinceSlipRateBelowThreshold[static_cast<size_t>(j * SimCount) + sim]),
              localMin[sim]);
        }
      }
    }
  }

  for (std::size_t sim = 0; sim < SimCount; ++sim) {
    minTimeSinceSlipRateBelowThreshold[sim] =
        std::min(localMin[sim], minTimeSinceSlipRateBelowThreshold[sim]);
  }

  for (std::size_t sim = 0; sim < SimCount; ++sim) {
    energies.energy(TotalFrictionalWork, sim) += totalFrictionalWork[sim];
    energies.energy(StaticFrictionalWork, sim) += staticFrictionalWork[sim];
    energies.energy(SeismicMoment, sim) += seismicMoment[sim];
    energies.energy(Potency, sim) += potency[sim];
  }
}

/// Adds the energies of the cells of `layer`, which compute in the configuration `Cfg`.
template <typename Cfg>
void addVolumeEnergies(const LTS::Layer& layer,
                       const std::vector<Element>& elements,
                       const std::vector<Vertex>& vertices,
                       double g,
                       bool isPlasticityEnabled,
                       EnergiesStorage& energies) {
  using real = Real<Cfg>;
  using MaterialT = model::MaterialOf<Cfg>;
  using EnergyComputeT = model::EnergyCompute<MaterialT>;
  constexpr auto Variant = configIdOf<Cfg>();
  constexpr auto SimCount = Cfg::NumSimulations;

  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadraturePointsTet = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const auto quadratureTet = seissol::quadrature::simplexRule<3>(QuadPolyDegree);
  const auto& quadratureWeightsTet = quadratureTet.second;

  // Note: Default(none) is not possible, clang requires data sharing attribute for g, gcc forbids
  // it
  const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
  const auto* cellInformationData = layer.var<LTS::CellInformation>();
  const auto* faceDisplacementsData = layer.var<LTS::FaceDisplacements>(Cfg());
  const auto* materialData = layer.var<LTS::MaterialData>(Cfg());
  const auto* boundaryMappingData = layer.var<LTS::BoundaryMapping>(Cfg());
  const auto* pstrainData = layer.var<LTS::PStrain>(Cfg());
  const auto* dofsData = layer.var<LTS::Dofs>(Cfg());
  const auto* energyData = layer.var<LTS::EnergyData>(Cfg());
  // only allocated for materials with anelastic variables
  const auto* dofsAneData = layer.var<LTS::DofsAne>(Cfg());

  constexpr auto EnergyCountSingle = EnergyComputeT::EnergyCount;
  constexpr auto EnergyCount = EnergyCountSingle * SimCount;

  double energyValues[EnergyCount]{};
  double localPlasticMoment[SimCount]{};
  double localGravitationalPotentialEnergy[SimCount]{};

#if !NVHPC_AVOID_OMP
#pragma omp parallel for schedule(static)                                                          \
    reduction(+ : localGravitationalPotentialEnergy[ : SimCount],                                  \
                  energyValues[ : EnergyCount],                                                    \
                  localPlasticMoment[ : SimCount])                                                 \
    shared(elements, vertices, quadratureWeightsTet)
#endif
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    if (secondaryInformation[cell].duplicate > 0) {
      // skip duplicate cells
      continue;
    }
    const auto elementId = secondaryInformation[cell].meshId;
    const double volume = MeshTools::volume(elements[elementId], vertices);

    // NOLINTNEXTLINE
    const auto& material = materialData[cell];
    const auto& cellInformation = cellInformationData[cell];
    const auto& faceDisplacements = faceDisplacementsData[cell];

    // Needed to weight the integral.
    const auto jacobiDet = 6 * volume;

    alignas(Alignment) real linData[tensor::momentQ<Cfg>::size()];
    auto lin = init::momentQ<Cfg>::view::create(linData);
    // cell integral of Q: momentQ(0, J) == \int_{T_ref} Q_J
    runtime::kernel::momentQCompute krnl;
    krnl.momentQ = runtime::init::momentQ::view(Variant, linData);
    krnl.Q = runtime::init::Q::view(Variant, dofsData[cell]);
    krnl.execute(Variant);

    alignas(Alignment) real quadData[tensor::momentQQ<Cfg>::size()];
    auto quad = init::momentQQ<Cfg>::view::create(quadData);
    // second moments of Q: momentQQ(I, J) == \int_{T_ref} Q_I Q_J
    runtime::kernel::momentQQCompute krnl2;
    krnl2.momentQQ = runtime::init::momentQQ::view(Variant, quadData);
    krnl2.Q = runtime::init::Q::view(Variant, dofsData[cell]);
    krnl2.execute(Variant);

    const auto moments = EnergyComputeT::template computeMoments<Cfg>(
        dofsData[cell], dofsAneData != nullptr ? dofsAneData[cell] : nullptr);

    for (size_t sim = 0; sim < SimCount; sim++) {

      auto linSub = multisim::simtensor<Cfg>(lin, sim);
      auto quadSub = multisim::simtensor<Cfg>(quad, sim);

      // assume _constant_ material over a cell (will need adjustments for e.g. #1297)

      const auto localValues = EnergyComputeT::template computeEnergies<Cfg>(
          material, energyData[cell], linSub, quadSub, moments, sim);

      for (std::size_t i = 0; i < localValues.size(); ++i) {
        energyValues[localValues.size() * sim + i] += jacobiDet * localValues[i];
      }
    }

    constexpr auto UIdx = MaterialT::VelocityOffset;

    const auto& boundaryMappings = boundaryMappingData[cell];
    // Compute the gravitational potential energy
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (cellInformation.faceTypes[face] != FaceType::FreeSurfaceGravity) {
        continue;
      }

      // Displacements are stored in face-aligned coordinate system.
      // We need to rotate it to the global coordinate system.
      const auto& boundaryMapping = boundaryMappings[face];
      auto tinv = init::Tinv<Cfg>::view::create(boundaryMapping.dataTinv);
      alignas(Alignment)
          real rotateDisplacementToFaceNormalData[init::displacementRotationMatrix<Cfg>::Size];

      auto rotateDisplacementToFaceNormal =
          init::displacementRotationMatrix<Cfg>::view::create(rotateDisplacementToFaceNormalData);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          rotateDisplacementToFaceNormal(i, j) = tinv(i + UIdx, j + UIdx);
        }
      }

      const auto* curFaceDisplacementsData = faceDisplacements[face];

      // See for example (Saito, Tsunami generation and propagation, 2019) section 3.2.3 for
      // derivation.
      //
      // The rotation into the global frame has to happen *before* squaring, hence the
      // two-step approach: the kernel produces the modal coefficients of the rotated
      // displacement, and the quadratic form against M2 is evaluated here.

      alignas(Alignment) std::array<real, tensor::faceDisplacementSquared<Cfg>::Size>
          faceDisplacementSquared{};
      {
        runtime::kernel::faceDisplacementSquaredCompute evalKrnl;
        evalKrnl.rotatedFaceDisplacement =
            runtime::init::rotatedFaceDisplacement::view(Variant, curFaceDisplacementsData);
        evalKrnl.faceDisplacementSquared =
            runtime::init::faceDisplacementSquared::view(Variant, faceDisplacementSquared.data());
        evalKrnl.displacementRotationMatrix = runtime::init::displacementRotationMatrix::view(
            Variant, rotateDisplacementToFaceNormalData);
        evalKrnl.execute(Variant);
      }

      const auto squaredViewFused =
          init::faceDisplacementSquared<Cfg>::view::create(faceDisplacementSquared.data());

      const auto surface = MeshTools::surface(elements[elementId], face, vertices);
      const auto rho = material.getDensity();

      for (size_t sim = 0; sim < SimCount; sim++) {
        const auto squaredView = multisim::simtensor<Cfg>(squaredViewFused, sim);

        // contains an elided 0.5 * 2.0 (1/2 due to energy; 2 due to surface)
        localGravitationalPotentialEnergy[sim] += rho * g * surface * squaredView(0);
      }
    }

    if (isPlasticityEnabled) {
      // plastic moment
      const real* pstrainCell = pstrainData[cell];
      const double mu = material.getMuBar();

      // integrating over all collocation points suffices
      const real* __restrict qEta = &pstrainCell[tensor::QStressNodal<Cfg>::size()];

      alignas(Alignment) real qEtaQuad[tensor::QEtaNodalProject<Cfg>::size()]{};

      runtime::kernel::plProject krnl;
      krnl.QEtaNodal = runtime::init::QEtaNodal::view(Variant, qEta);
      krnl.QEtaNodalProject = runtime::init::QEtaNodalProject::view(Variant, qEtaQuad);
      krnl.execute(Variant);

      // go through the view: QEtaNodalProject is padded (at order 6, its 343 points take up 344
      // entries), and for fused simulations the simulation index leads and may be padded as well
      static_assert(tensor::QEtaNodalProject<Cfg>::Shape[multisim::BasisDim<Cfg>] ==
                    NumQuadraturePointsTet);
      auto qEtaQuadView = init::QEtaNodalProject<Cfg>::view::create(qEtaQuad);
      for (size_t sim = 0; sim < SimCount; ++sim) {
        const auto qEtaQuadSim = multisim::simtensor<Cfg>(qEtaQuadView, sim);
        double pMoment = 0;
        for (size_t qp = 0; qp < NumQuadraturePointsTet; ++qp) {
          pMoment += quadratureWeightsTet[qp] * qEtaQuadSim(qp);
        }
        localPlasticMoment[sim] += mu * jacobiDet * pMoment;
      }
    }
  }

  for (std::size_t sim = 0; sim < SimCount; ++sim) {
    for (std::size_t i = 0; i < EnergyComputeT::EnergyCount; ++i) {
      const auto& descriptor = EnergyComputeT::Energies[i];
      energies.energy(descriptor.name, sim) += energyValues[sim * EnergyComputeT::EnergyCount + i];
    }

    energies.energy(PlasticMoment, sim) += localPlasticMoment[sim];
    energies.energy(GravitationalPotentialEnergy, sim) += localGravitationalPotentialEnergy[sim];
  }
}

} // namespace

const SIUnit& siUnit(model::EnergyUnit unit) {
  switch (unit) {
  case model::EnergyUnit::Energy:
    return UnitEnergy;
  case model::EnergyUnit::Power:
    return UnitPower;
  case model::EnergyUnit::Moment:
    return UnitMoment;
  case model::EnergyUnit::Momentum:
    return UnitMomentum;
  case model::EnergyUnit::Scalar:
    return UnitScalar;
  }
  logError() << "Unhandled energy unit.";
  return UnitScalar;
}

void EnergiesStorage::setSimcount(size_t count) { simcount_ = count; }

size_t EnergiesStorage::addEnergy(const model::EnergyDescriptor& descriptor) {
  if (descriptor.name.empty()) {
    logError() << "Attempted to register an energy without a name. This usually means that an"
               << "EnergyCompute specialization declares more energies than it names.";
  }
  if (handles_.find(descriptor.name) != handles_.end()) {
    logError() << "Energy" << std::string(descriptor.name) << "registered twice.";
  }
  const auto index = descriptors_.size();
  handles_.emplace(std::string(descriptor.name), index);
  descriptors_.emplace_back(descriptor);
  values_.resize(values_.size() + simcount_, 0);
  return index;
}

double& EnergiesStorage::energy(size_t handle, size_t sim) {
  return values_[simcount_ * handle + sim];
}

[[nodiscard]] double EnergiesStorage::energy(size_t handle, size_t sim) const {
  return values_[simcount_ * handle + sim];
}

[[nodiscard]] size_t EnergiesStorage::handleOf(std::string_view name) const {
  const auto it = handles_.find(name);
  if (it == handles_.end()) {
    // Deliberately fatal: silently returning zero turns a typo into a
    // plausible-looking number that nobody notices.
    logError() << "Unknown energy" << std::string(name).c_str()
               << "-- use EnergiesStorage::has() to test for optional quantities.";
  }
  return it->second;
}

double& EnergiesStorage::energy(std::string_view name, size_t sim) {
  return energy(handleOf(name), sim);
}

[[nodiscard]] double EnergiesStorage::energy(std::string_view name, size_t sim) const {
  return energy(handleOf(name), sim);
}

[[nodiscard]] bool EnergiesStorage::has(std::string_view name) const {
  return handles_.find(name) != handles_.end();
}

[[nodiscard]] const std::vector<model::EnergyDescriptor>& EnergiesStorage::descriptors() const {
  return descriptors_;
}

std::vector<double>& EnergiesStorage::values() { return values_; }

void EnergiesStorage::reset() { std::fill(values_.begin(), values_.end(), 0); }

void EnergyOutput::init(
    const DynamicRupture::Storage& newDynRuptTree,
    const seissol::geometry::MeshReader& newMeshReader,
    const LTS::Storage& newStorage,
    bool newIsPlasticityEnabled,
    const std::string& outputFileNamePrefix,
    const seissol::initializer::parameters::EnergyOutputParameters& parameters) {
  if (parameters.enabled && parameters.interval > 0) {
    isEnabled_ = true;
  } else {
    return;
  }
  const auto rank = Mpi::mpi.rank();
  logInfo() << "Initializing energy output.";

  energyOutputInterval_ = parameters.interval;
  isFileOutputEnabled_ = rank == 0;
  isTerminalOutputEnabled_ = parameters.terminalOutput && (rank == 0);
  terminatorMaxTimePostRupture_ = parameters.terminatorMaxTimePostRupture;
  terminatorMomentRateThreshold_ = parameters.terminatorMomentRateThreshold;
  // The slip-rate terminator is active exactly when a finite post-rupture time was
  // configured. The comparison used to be the other way round, which enabled the
  // check precisely when the user had switched the terminator off.
  isCheckAbortCriteraSlipRateEnabled_ =
      (terminatorMaxTimePostRupture_ < std::numeric_limits<double>::max());
  isCheckAbortCriteraMomentRateEnabled_ = (terminatorMomentRateThreshold_ > 0);
  computeVolumeEnergiesEveryOutput_ = parameters.computeVolumeEnergiesEveryOutput;
  outputFileName_ = outputFileNamePrefix + "-energy.csv";

  drStorage_ = &newDynRuptTree;
  meshReader_ = &newMeshReader;
  ltsStorage_ = &newStorage;

  isPlasticityEnabled_ = newIsPlasticityEnabled;

  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  setSyncInterval(parameters.interval);

  // Every rank registers the same energies, in the same order: the ones of the materials of the
  // configurations the run has, each once, also when several of these configurations share a
  // material.
  const auto configs = ltsStorage_->configs();
  simulationCount_ = 1;
  for (const auto config : configs) {
    dispatchConfig(config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      simulationCount_ = std::max<std::size_t>(simulationCount_, Cfg::NumSimulations);
    });
  }

  energiesStorage_.setSimcount(simulationCount_);
  minTimeSinceSlipRateBelowThreshold_.assign(simulationCount_, 0.0);
  minTimeSinceMomentRateBelowThreshold_.assign(simulationCount_, 0.0);
  seismicMomentPrevious_.assign(simulationCount_, 0.0);

  for (const auto& descriptor : GlobalEnergies) {
    energiesStorage_.addEnergy(descriptor);
  }

  for (const auto config : configs) {
    dispatchConfig(config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      for (const auto& descriptor : model::EnergyCompute<model::MaterialOf<Cfg>>::Energies) {
        if (!energiesStorage_.has(descriptor.name)) {
          energiesStorage_.addEnergy(descriptor);
        }
      }
    });
  }
}

void EnergyOutput::syncPoint(double time) {
  assert(isEnabled_);
  const auto rank = Mpi::mpi.rank();
  logInfo() << "Writing energy output at time" << time;

  seissolInstance_.dofSync().syncDofs(time);

  computeEnergies();
  reduceEnergies();
  if (isCheckAbortCriteraSlipRateEnabled_) {
    reduceMinTimeSinceSlipRateBelowThreshold();
  }
  if ((rank == 0) && isCheckAbortCriteraMomentRateEnabled_) {
    for (size_t sim = 0; sim < simulationCount_; sim++) {
      const double seismicMomentRate =
          (energiesStorage_.energy(SeismicMoment, sim) - seismicMomentPrevious_[sim]) /
          energyOutputInterval_;
      seismicMomentPrevious_[sim] = energiesStorage_.energy(SeismicMoment, sim);
      if (time > 0 && seismicMomentRate < terminatorMomentRateThreshold_) {
        minTimeSinceMomentRateBelowThreshold_[sim] += energyOutputInterval_;
      } else {
        minTimeSinceMomentRateBelowThreshold_[sim] = 0.0;
      }
    }
  }
  if (isTerminalOutputEnabled_) {
    printEnergies();
  }
  if (isCheckAbortCriteraSlipRateEnabled_) {
    checkAbortCriterion(minTimeSinceSlipRateBelowThreshold_, "All slip rates are");
  }
  if (isCheckAbortCriteraMomentRateEnabled_) {
    checkAbortCriterion(minTimeSinceMomentRateBelowThreshold_, "The seismic moment rate is");
  }

  if (isFileOutputEnabled_) {
    writeEnergies(time);
  }
  ++outputId_;
  logInfo() << "Writing energy output at time" << time << "Done.";
}

void EnergyOutput::simulationStart(std::optional<double> checkpointTime) {
  if (isFileOutputEnabled_) {
    // a run resuming from a checkpoint keeps the energies up to it, as the other outputs do
    if (checkpointTime.has_value()) {
      io::writer::file::backUpFile(outputFileName_);
    }
    std::size_t nameWidth = 0;
    for (const auto& descriptor : energiesStorage_.descriptors()) {
      nameWidth = std::max(nameWidth, descriptor.name.size());
    }
    table_.emplace("energy");
    table_->addColumn<double>("time");
    table_->addTextColumn("variable", nameWidth);
    table_->addColumn<std::uint64_t>("simulation_index");
    table_->addColumn<double>("measurement");
  }
  syncPoint(checkpointTime.value_or(0));
}

EnergyOutput::~EnergyOutput() = default;

void EnergyOutput::computeDynamicRuptureEnergies() {
  std::fill(minTimeSinceSlipRateBelowThreshold_.begin(),
            minTimeSinceSlipRateBelowThreshold_.end(),
            std::numeric_limits<double>::max());

  auto& memoryManager = seissolInstance_.memoryManager();
  for (const auto& layer : drStorage_->leaves()) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      addFaultEnergies<Cfg>(layer,
                            *memoryManager.globalData<Cfg>().onHost,
                            energiesStorage_,
                            minTimeSinceSlipRateBelowThreshold_);
    });
  }
}

void EnergyOutput::computeVolumeEnergies() {
  const std::vector<Element>& elements = meshReader_->getElements();
  const std::vector<Vertex>& vertices = meshReader_->getVertices();

  const auto g = seissolInstance_.gravitationSetup().acceleration;

  for (const auto& layer : ltsStorage_->leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      addVolumeEnergies<Cfg>(layer, elements, vertices, g, isPlasticityEnabled_, energiesStorage_);
    });
  }
}

void EnergyOutput::computeEnergies() {
  energiesStorage_.reset();
  if (shouldComputeVolumeEnergies()) {
    computeVolumeEnergies();
  }
  computeDynamicRuptureEnergies();
}

void EnergyOutput::reduceEnergies() {
  const auto& comm = Mpi::mpi.comm();
  MPI_Allreduce(MPI_IN_PLACE,
                energiesStorage_.values().data(),
                static_cast<int>(energiesStorage_.values().size()),
                MPI_DOUBLE,
                MPI_SUM,
                comm);
}

void EnergyOutput::reduceMinTimeSinceSlipRateBelowThreshold() {
  const auto& comm = Mpi::mpi.comm();
  MPI_Allreduce(MPI_IN_PLACE,
                minTimeSinceSlipRateBelowThreshold_.data(),
                static_cast<int>(minTimeSinceSlipRateBelowThreshold_.size()),
                Mpi::castToMpiType<double>(),
                MPI_MIN,
                comm);
}

void EnergyOutput::printEnergies() {
  const auto outputPrecision =
      seissolInstance_.parameters().output.energyParameters.terminalPrecision;

  std::vector<std::pair<std::size_t, std::string>> infnan;

  const auto shouldPrint = [](double thresholdValue) { return std::abs(thresholdValue) > 1.e-20; };
  for (size_t sim = 0; sim < simulationCount_; sim++) {
    const std::string fusedPrefix = simulationCount_ > 1 ? "[" + std::to_string(sim) + "]" : "";

    const auto printValue = [&](double value, const SIUnit& unit) {
      return unit.formatScientific(value, {}, outputPrecision);
    };
    const auto magnitude = [&](double moment) { return 2.0 / 3.0 * std::log10(moment) - 6.07; };

    // Energies sharing a group are summed and reported on one line, driven by the
    // member that carries the heading. A single-member group prints just the
    // total; a larger one appends each member's share. This is what lets the
    // viscoelastic branch energy join the elastic group -- reporting the
    // kinetic/potential split without it would understate the potential share.
    const auto printGroup = [&](const model::EnergyDescriptor& labelled) {
      const auto& descriptors = energiesStorage_.descriptors();
      double total = 0.0;
      for (const auto& member : descriptors) {
        if (member.group == labelled.group) {
          total += energiesStorage_.energy(member.name, sim);
        }
      }
      // guard before dividing, so an empty field does not produce 0/0
      if (!shouldPrint(total)) {
        return;
      }

      std::ostringstream shares;
      shares << std::setprecision(outputPrecision);
      for (const auto& member : descriptors) {
        if (member.group != labelled.group || member.shortLabel.empty()) {
          continue;
        }
        shares << " , " << member.shortLabel << " "
               << (energiesStorage_.energy(member.name, sim) / total * 100.0) << " %";
      }

      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << std::string(labelled.groupLabel).c_str()
                << printValue(total, siUnit(labelled.unit)).c_str() << shares.str().c_str();
    };

    const auto seismicMoment = energiesStorage_.energy(SeismicMoment, sim);
    if (shouldComputeVolumeEnergies()) {
      // Every group is printed generically, in descriptor registration order.
      // Anything needing more than a total and per-member shares -- the plastic
      // moment with its equivalent magnitude, the momentum triple, frictional
      // work -- gets a bespoke block below.
      for (const auto& descriptor : energiesStorage_.descriptors()) {
        if (!descriptor.groupLabel.empty()) {
          printGroup(descriptor);
        }
      }

      const auto plasticMoment = energiesStorage_.energy(PlasticMoment, sim);
      if (shouldPrint(plasticMoment)) {
        const auto ratioPlasticMoment = 100.0 * plasticMoment / (plasticMoment + seismicMoment);
        logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                  << "Plastic moment:" << printValue(plasticMoment, UnitMoment).c_str()
                  << ", equivalent Mw:" << magnitude(plasticMoment)
                  << ", of total moment:" << ratioPlasticMoment << "%";
      }

      if (std::all_of(MomentumComponents.begin(),
                      MomentumComponents.end(),
                      [&](std::string_view name) { return energiesStorage_.has(name); })) {
        logInfo()
            << std::setprecision(outputPrecision) << fusedPrefix.c_str() << " Total momentum: X"
            << printValue(energiesStorage_.energy(MomentumComponents[0], sim), UnitMomentum).c_str()
            << ", Y"
            << printValue(energiesStorage_.energy(MomentumComponents[1], sim), UnitMomentum).c_str()
            << ", Z"
            << printValue(energiesStorage_.energy(MomentumComponents[2], sim), UnitMomentum)
                   .c_str();
      }
    } else {
      logInfo() << "Volume energies skipped at this step";
    }

    const auto totalFrictionalWork = energiesStorage_.energy(TotalFrictionalWork, sim);
    if (shouldPrint(totalFrictionalWork)) {
      const auto staticFrictionalWork = energiesStorage_.energy(StaticFrictionalWork, sim);
      const auto ratio1 = staticFrictionalWork / totalFrictionalWork * 100.0;
      const auto ratio2 =
          (totalFrictionalWork - staticFrictionalWork) / totalFrictionalWork * 100.0;
      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << "Frictional work: " << printValue(totalFrictionalWork, UnitEnergy).c_str()
                << ", static" << ratio1 << "% , radiated" << ratio2 << "%";

      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << "Seismic moment (without plasticity):"
                << printValue(seismicMoment, UnitMoment).c_str()
                << ", Mw:" << magnitude(seismicMoment);
    }

    for (const auto& descriptor : energiesStorage_.descriptors()) {
      if (!std::isfinite(energiesStorage_.energy(descriptor.name, sim))) {
        infnan.emplace_back(sim, std::string(descriptor.name));
      }
    }
  }

  if (!infnan.empty()) {
    logError() << "Detected Inf/NaN in energies. Aborting. Inf/NaN values:" << infnan;
  }
}

void EnergyOutput::checkAbortCriterion(const std::vector<double>& timeSinceThreshold,
                                       const std::string& prefixMessage) {
  // A simulation counts as "ready to abort" once it has been below the threshold
  // for longer than the configured time. Simulations that have not reached the
  // threshold at all are still running and must block the abort; simulations that
  // never entered the observation window (timeSinceThreshold == 0 or infinite)
  // carry no information and must *not* block it, or a single such simulation
  // would keep a fused run alive indefinitely.
  size_t abortCount = 0;
  size_t decidableCount = 0;
  for (size_t sim = 0; sim < simulationCount_; sim++) {
    if ((timeSinceThreshold[sim] > 0) and
        (timeSinceThreshold[sim] < std::numeric_limits<double>::infinity())) {
      ++decidableCount;
      if (static_cast<double>(timeSinceThreshold[sim]) < terminatorMaxTimePostRupture_) {
        logInfo() << prefixMessage.c_str() << "below threshold since" << timeSinceThreshold[sim]
                  << "s; in simulation: " << sim
                  << "(lower than the abort criteria: " << terminatorMaxTimePostRupture_ << "s)";
      } else {
        logInfo() << prefixMessage.c_str() << "below threshold since" << timeSinceThreshold[sim]
                  << "s; in simulation: " << sim
                  << "(greater than the abort criteria: " << terminatorMaxTimePostRupture_ << "s)";
        ++abortCount;
      }
    }
  }

  bool abort = (decidableCount > 0) and (abortCount == decidableCount);
  const auto& comm = Mpi::mpi.comm();
  MPI_Bcast(reinterpret_cast<void*>(&abort), 1, MPI_CXX_BOOL, 0, comm);
  if (abort) {
    seissolInstance_.simulator().abort();
  }
}

void EnergyOutput::writeEnergies(double time) {
  // iterate the descriptors, not the name->handle map: the map is ordered
  // alphabetically, the descriptors in registration order
  const auto& descriptors = energiesStorage_.descriptors();
  for (std::size_t handle = 0; handle < descriptors.size(); ++handle) {
    for (size_t sim = 0; sim < simulationCount_; sim++) {
      table_->addCell<double>(time);
      table_->addText(std::string(descriptors[handle].name));
      table_->addCell<std::uint64_t>(sim);
      table_->addCell<double>(energiesStorage_.energy(handle, sim));
    }
  }
  table_->appendFile(outputFileName_);
}

bool EnergyOutput::shouldComputeVolumeEnergies() const {
  return outputId_ % computeVolumeEnergiesEveryOutput_ == 0;
}

} // namespace seissol::writer
