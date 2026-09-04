// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EnergyOutput.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "Equations/Energy.h"
#include "Equations/EnergyBase.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/PreProcessorMacros.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Modules/Modules.h"
#include "Monitoring/Unit.h"
#include "Numerical/Quadrature.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <ios>
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
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

GENERATE_HAS_MEMBER(vInv)
GENERATE_HAS_MEMBER(evalAtQP)

namespace seissol::writer {

namespace {

std::array<real, multisim::NumSimulations>
    computeStaticWork(const real* degreesOfFreedomPlus,
                      const real* degreesOfFreedomMinus,
                      const DRFaceInformation& faceInfo,
                      const DRGodunovData& godunovData,
                      const real slip[seissol::tensor::slipInterpolated::size()],
                      const GlobalData* global) {
  real points[seissol::kernels::NumSpaceQuadraturePoints][2];
  alignas(Alignment) real spaceWeights[seissol::kernels::NumSpaceQuadraturePoints];
  seissol::quadrature::TriangleQuadrature(points, spaceWeights, ConvergenceOrder + 1);

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints krnl;
  krnl.V3mTo2n = global->faceToNodalMatrices;

  alignas(PagesizeStack) real qInterpolatedPlus[tensor::QInterpolatedPlus::size()];
  alignas(PagesizeStack) real qInterpolatedMinus[tensor::QInterpolatedMinus::size()];
  alignas(Alignment) real tractionInterpolated[tensor::tractionInterpolated::size()];
  alignas(Alignment) real qPlus[tensor::Q::size()];
  alignas(Alignment) real qMinus[tensor::Q::size()];

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

  dynamicRupture::kernel::computeTractionInterpolated trKrnl;
  trKrnl.tractionPlusMatrix = godunovData.tractionPlusMatrix;
  trKrnl.tractionMinusMatrix = godunovData.tractionMinusMatrix;
  trKrnl.QInterpolatedPlus = qInterpolatedPlus;
  trKrnl.QInterpolatedMinus = qInterpolatedMinus;
  trKrnl.tractionInterpolated = tractionInterpolated;
  trKrnl.execute();

  alignas(Alignment) real staticFrictionalWork[tensor::staticFrictionalWork::size()]{};

  dynamicRupture::kernel::accumulateStaticFrictionalWork feKrnl;
  feKrnl.slipInterpolated = slip;
  feKrnl.tractionInterpolated = tractionInterpolated;
  feKrnl.spaceWeights = spaceWeights;
  feKrnl.staticFrictionalWork = staticFrictionalWork;
  feKrnl.minusSurfaceArea = -0.5 * godunovData.doubledSurfaceArea;
  feKrnl.execute();

  std::array<real, multisim::NumSimulations> frictionalWorkReturn{};
  std::copy_n(staticFrictionalWork, multisim::NumSimulations, frictionalWorkReturn.begin());
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
constexpr std::string_view GravitationalEnergy = "gravitational_energy";
constexpr std::string_view SeismicMoment = "seismic_moment";
constexpr std::string_view TotalFrictionalWork = "total_frictional_work";
constexpr std::string_view StaticFrictionalWork = "static_frictional_work";
constexpr std::string_view Potency = "potency";

constexpr std::array GlobalEnergies{
    model::EnergyDescriptor{PlasticMoment, model::EnergyUnit::Moment, {}, {}, {}},
    model::EnergyDescriptor{GravitationalEnergy,
                            model::EnergyUnit::Energy,
                            "gravitational",
                            "Gravitational energy:",
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
    GlobalData* newGlobal,
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

  global_ = newGlobal;
  drStorage_ = &newDynRuptTree;
  meshReader_ = &newMeshReader;
  ltsStorage_ = &newStorage;

  isPlasticityEnabled_ = newIsPlasticityEnabled;

  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  setSyncInterval(parameters.interval);

  energiesStorage_.setSimcount(multisim::NumSimulations);

  for (const auto& descriptor : GlobalEnergies) {
    energiesStorage_.addEnergy(descriptor);
  }

  for (const auto& descriptor : model::EnergyCompute<model::MaterialT>::Energies) {
    energiesStorage_.addEnergy(descriptor);
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
    for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
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
    out_.open(outputFileName_);
    out_ << std::scientific;
    out_ << std::setprecision(std::numeric_limits<double>::max_digits10);
    writeHeader();
  }
  syncPoint(checkpointTime.value_or(0));
}

EnergyOutput::~EnergyOutput() = default;

void EnergyOutput::computeDynamicRuptureEnergies() {
  constexpr auto SimCount = multisim::NumSimulations;

  double totalFrictionalWork[SimCount]{};
  double staticFrictionalWork[SimCount]{};
  double seismicMoment[SimCount]{};
  double potency[SimCount]{};

  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    minTimeSinceSlipRateBelowThreshold_[sim] = std::numeric_limits<double>::max();
  }

  for (const auto& layer : drStorage_->leaves()) {

    real* const* timeDofsPlus = layer.var<DynamicRupture::TimeDerivativePlus>();
    real* const* timeDofsMinus = layer.var<DynamicRupture::TimeDerivativeMinus>();

    const auto* godunovData = layer.var<DynamicRupture::GodunovData>();
    const auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
    const auto* drEnergyOutput = layer.var<DynamicRupture::DREnergyOutputVar>();
    const auto* waveSpeedsPlus = layer.var<DynamicRupture::WaveSpeedsPlus>();
    const auto* waveSpeedsMinus = layer.var<DynamicRupture::WaveSpeedsMinus>();
    const auto layerSize = layer.size();

#if !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(+ : totalFrictionalWork[ : SimCount],                           \
                                       staticFrictionalWork[ : SimCount],                          \
                                       seismicMoment[ : SimCount],                                 \
                                       potency[ : SimCount]) default(none)                         \
    shared(layerSize,                                                                              \
               drEnergyOutput,                                                                     \
               faceInformation,                                                                    \
               timeDofsMinus,                                                                      \
               timeDofsPlus,                                                                       \
               godunovData,                                                                        \
               waveSpeedsPlus,                                                                     \
               waveSpeedsMinus,                                                                    \
               SimCount)
#endif
    for (std::size_t i = 0; i < layerSize; ++i) {
      if (faceInformation[i].plusSideOnThisRank) {
        const auto staticFrictionalWorkIncrease = computeStaticWork(timeDofsPlus[i],
                                                                    timeDofsMinus[i],
                                                                    faceInformation[i],
                                                                    godunovData[i],
                                                                    drEnergyOutput[i].slip,
                                                                    global_);

        const double muPlus = waveSpeedsPlus[i].density * waveSpeedsPlus[i].sWaveVelocity *
                              waveSpeedsPlus[i].sWaveVelocity;
        const double muMinus = waveSpeedsMinus[i].density * waveSpeedsMinus[i].sWaveVelocity *
                               waveSpeedsMinus[i].sWaveVelocity;
        const double mu = 2.0 * muPlus * muMinus / (muPlus + muMinus);

#pragma omp simd
        for (size_t sim = 0; sim < SimCount; sim++) {
          staticFrictionalWork[sim] += staticFrictionalWorkIncrease[sim];
          for (std::size_t j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints; ++j) {
            totalFrictionalWork[sim] += drEnergyOutput[i].frictionalEnergy[j * SimCount + sim];
          }
          double potencyIncrease = 0.0;
          for (std::size_t k = 0; k < seissol::dr::misc::NumBoundaryGaussPoints; ++k) {
            potencyIncrease += drEnergyOutput[i].accumulatedSlip[k * SimCount + sim] *
                               init::quadweights::Values[k];
          }
          potencyIncrease *= godunovData[i].doubledSurfaceArea;
          potency[sim] += potencyIncrease;
          seismicMoment[sim] += potencyIncrease * mu;
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
          for (std::size_t j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints; ++j) {
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
      minTimeSinceSlipRateBelowThreshold_[sim] =
          std::min(localMin[sim], minTimeSinceSlipRateBelowThreshold_[sim]);
    }
  }

  for (std::size_t sim = 0; sim < SimCount; ++sim) {
    energiesStorage_.energy(TotalFrictionalWork, sim) += totalFrictionalWork[sim];
    energiesStorage_.energy(StaticFrictionalWork, sim) += staticFrictionalWork[sim];
    energiesStorage_.energy(SeismicMoment, sim) += seismicMoment[sim];
    energiesStorage_.energy(Potency, sim) += potency[sim];
  }
}

void EnergyOutput::computeVolumeEnergies() {
  const std::vector<Element>& elements = meshReader_->getElements();
  const std::vector<Vertex>& vertices = meshReader_->getVertices();

  const auto g = seissolInstance_.gravitationSetup().acceleration;

  constexpr auto QuadPolyDegree = ConvergenceOrder + 1;
  constexpr auto NumQuadraturePointsTet = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  double quadraturePointsTet[NumQuadraturePointsTet][3]{};
  double quadratureWeightsTet[NumQuadraturePointsTet]{};
  seissol::quadrature::TetrahedronQuadrature(
      quadraturePointsTet, quadratureWeightsTet, QuadPolyDegree);

  // Note: Default(none) is not possible, clang requires data sharing attribute for g, gcc forbids
  // it
  for (const auto& layer : ltsStorage_->leaves(Ghost)) {
    const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    const auto* cellInformationData = layer.var<LTS::CellInformation>();
    const auto* faceDisplacementsData = layer.var<LTS::FaceDisplacements>();
    const auto* materialData = layer.var<LTS::MaterialData>();
    const auto* boundaryMappingData = layer.var<LTS::BoundaryMapping>();
    const auto* pstrainData = layer.var<LTS::PStrain>();
    const auto* dofsData = layer.var<LTS::Dofs>();
    const auto* energyData = layer.var<LTS::EnergyData>();
    // only allocated for materials with anelastic variables
    const auto* dofsAneData = layer.var<LTS::DofsAne>();

    constexpr auto SimCount = multisim::NumSimulations;
    constexpr auto EnergyCountSingle = model::EnergyCompute<model::MaterialT>::EnergyCount;
    constexpr auto EnergyCount = EnergyCountSingle * SimCount;

    double energyValues[EnergyCount]{};
    double localPlasticMoment[SimCount]{};
    double localGravitationalEnergy[SimCount]{};

#if !NVHPC_AVOID_OMP
#pragma omp parallel for schedule(static) reduction(+ : localGravitationalEnergy[ : SimCount],     \
                                                        energyValues[ : EnergyCount],              \
                                                        localPlasticMoment[ : SimCount])           \
    shared(elements, vertices, global_)
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

      alignas(Alignment) real linData[tensor::momentQ::size()];
      auto lin = init::momentQ::view::create(linData);
      // cell integral of Q: momentQ(0, J) == \int_{T_ref} Q_J
      kernel::momentQCompute krnl;
      krnl.M3 = init::M3::Values;
      krnl.momentQ = linData;
      krnl.Q = dofsData[cell];
      krnl.execute();

      alignas(Alignment) real quadData[tensor::momentQQ::size()];
      auto quad = init::momentQQ::view::create(quadData);
      // second moments of Q: momentQQ(I, J) == \int_{T_ref} Q_I Q_J
      kernel::momentQQCompute krnl2;
      krnl2.M3 = init::M3::Values;
      krnl2.momentQQ = quadData;
      krnl2.Q = dofsData[cell];
      krnl2.execute();

      const auto moments = model::EnergyCompute<model::MaterialT>::computeMoments(
          dofsData[cell], dofsAneData != nullptr ? dofsAneData[cell] : nullptr);

      for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {

        auto linSub = multisim::simtensor(lin, sim);
        auto quadSub = multisim::simtensor(quad, sim);

        // assume _constant_ material over a cell (will need adjustments for e.g. #1297)

        const auto localValues = model::EnergyCompute<model::MaterialT>::computeEnergies(
            material, energyData[cell], linSub, quadSub, moments, sim);

        for (std::size_t i = 0; i < localValues.size(); ++i) {
          energyValues[localValues.size() * sim + i] += jacobiDet * localValues[i];
        }
      }

      constexpr auto UIdx = model::MaterialT::VelocityOffset;

      const auto& boundaryMappings = boundaryMappingData[cell];
      // Compute gravitational energy
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (cellInformation.faceTypes[face] != FaceType::FreeSurfaceGravity) {
          continue;
        }

        // Displacements are stored in face-aligned coordinate system.
        // We need to rotate it to the global coordinate system.
        const auto& boundaryMapping = boundaryMappings[face];
        auto tinv = init::Tinv::view::create(boundaryMapping.dataTinv);
        alignas(Alignment)
            real rotateDisplacementToFaceNormalData[init::displacementRotationMatrix::Size];

        auto rotateDisplacementToFaceNormal =
            init::displacementRotationMatrix::view::create(rotateDisplacementToFaceNormalData);
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

        alignas(Alignment) std::array<real, tensor::faceDisplacementSquared::Size>
            faceDisplacementSquared{};
        {
          seissol::kernel::faceDisplacementSquaredCompute evalKrnl;
          evalKrnl.rotatedFaceDisplacement = curFaceDisplacementsData;
          evalKrnl.M2 = init::M2::Values;
          evalKrnl.MV2nTo2m = nodal::init::MV2nTo2m::Values;
          evalKrnl.faceDisplacementSquared = faceDisplacementSquared.data();
          evalKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
          evalKrnl.execute();
        }

        const auto squaredViewFused =
            init::faceDisplacementSquared::view::create(faceDisplacementSquared.data());

        const auto surface = MeshTools::surface(elements[elementId], face, vertices);
        const auto rho = material.getDensity();

        for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
          const auto squaredView = multisim::simtensor(squaredViewFused, sim);

          // contains an elided 0.5 * 2.0 (1/2 due to energy; 2 due to surface)
          localGravitationalEnergy[sim] += rho * g * surface * squaredView(0);
        }
      }

      if (isPlasticityEnabled_) {
        // plastic moment
        const real* pstrainCell = pstrainData[cell];
        const double mu = material.getMuBar();

        // integrating over all collocation points suffices
        const real* __restrict qEta = &pstrainCell[tensor::QStressNodal::size()];

        alignas(Alignment) real qEtaQuad[tensor::QEtaNodalProject::size()]{};

        kernel::plProject krnl;
        set_evalAtQP(krnl, global_->evalAtQPMatrix);
        set_vInv(krnl, global_->vandermondeMatrixInverse);
        krnl.QEtaNodal = qEta;
        krnl.QEtaNodalProject = qEtaQuad;
        krnl.execute();

        // C-style array due to OpenMP
        double pMoment[multisim::NumSimulations]{};

#pragma omp simd reduction(+ : pMoment[ : multisim::NumSimulations])
        for (size_t qp = 0; qp < tensor::QEtaNodalProject::size(); ++qp) {
          pMoment[qp % multisim::NumSimulations] +=
              quadratureWeightsTet[qp / multisim::NumSimulations] * qEtaQuad[qp];
        }

        for (size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
          localPlasticMoment[sim] += mu * jacobiDet * pMoment[sim];
        }
      }
    }

    for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
      for (std::size_t i = 0; i < model::EnergyCompute<model::MaterialT>::EnergyCount; ++i) {
        const auto& descriptor = model::EnergyCompute<model::MaterialT>::Energies[i];
        energiesStorage_.energy(descriptor.name, sim) +=
            energyValues[sim * model::EnergyCompute<model::MaterialT>::EnergyCount + i];
      }

      energiesStorage_.energy(PlasticMoment, sim) += localPlasticMoment[sim];
      energiesStorage_.energy(GravitationalEnergy, sim) += localGravitationalEnergy[sim];
    }
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
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    const std::string fusedPrefix =
        multisim::MultisimEnabled ? "[" + std::to_string(sim) + "]" : "";

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

void EnergyOutput::checkAbortCriterion(
    const std::array<double, multisim::NumSimulations>& timeSinceThreshold,
    const std::string& prefixMessage) {
  // A simulation counts as "ready to abort" once it has been below the threshold
  // for longer than the configured time. Simulations that have not reached the
  // threshold at all are still running and must block the abort; simulations that
  // never entered the observation window (timeSinceThreshold == 0 or infinite)
  // carry no information and must *not* block it, or a single such simulation
  // would keep a fused run alive indefinitely.
  size_t abortCount = 0;
  size_t decidableCount = 0;
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
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

void EnergyOutput::writeHeader() {
  out_ << "time,variable,simulation_index,measurement" << std::endl;
}

void EnergyOutput::writeEnergies(double time) {
  // iterate the descriptors, not the name->handle map: the map is ordered
  // alphabetically, the descriptors in registration order
  const auto& descriptors = energiesStorage_.descriptors();
  for (std::size_t handle = 0; handle < descriptors.size(); ++handle) {
    for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
      out_ << time << "," << descriptors[handle].name << "," << sim << ","
           << energiesStorage_.energy(handle, sim) << '\n';
    }
  }
  out_.flush();
}

bool EnergyOutput::shouldComputeVolumeEnergies() const {
  return outputId_ % computeVolumeEnergiesEveryOutput_ == 0;
}

} // namespace seissol::writer
