// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EnergyOutput.h"

#include "DynamicRupture/Misc.h"
#include "Numerical/Quadrature.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include <Alignment.h>
#include <Common/ConfigHelper.h>
#include <Common/Constants.h>
#include <Equations/Datastructures.h>
#include <GeneratedCode/init.h>
#include <GeneratedCode/kernel.h>
#include <GeneratedCode/tensor.h>
#include <Geometry/MeshDefinition.h>
#include <Geometry/MeshTools.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/CellLocalInformation.h>
#include <Initializer/Parameters/OutputParameters.h>
#include <Initializer/PreProcessorMacros.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Model/CommonDatastructures.h>
#include <Modules/Modules.h>
#include <Solver/MultipleSimulations.h>
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
#include <mpi.h>
#include <optional>
#include <ostream>
#include <string>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#endif

namespace {
template <typename MaterialT>
constexpr bool VolumeEnergyApproximation = MaterialT::Type != model::MaterialType::Elastic &&
                                           MaterialT::Type != model::MaterialType::Viscoelastic &&
                                           MaterialT::Type != model::MaterialType::Acoustic;

template <typename Cfg>
std::array<Real<Cfg>, multisim::NumSimulations>
    computeStaticWork(const Real<Cfg>* degreesOfFreedomPlus,
                      const Real<Cfg>* degreesOfFreedomMinus,
                      const DRFaceInformation& faceInfo,
                      const DRGodunovData<Cfg>& godunovData,
                      const Real<Cfg>* slip,
                      const GlobalData* global) {
  // NOLINTNEXTLINE
  using real = Real<Cfg>;

  real points[seissol::kernels::NumSpaceQuadraturePoints<Cfg>][2];
  alignas(Alignment) real spaceWeights[seissol::kernels::NumSpaceQuadraturePoints<Cfg>];
  seissol::quadrature::TriangleQuadrature(points, spaceWeights, Cfg::ConvergenceOrder + 1);

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg> krnl;
  krnl.V3mTo2n = global->get<Cfg>().faceToNodalMatrices;

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

  dynamicRupture::kernel::computeTractionInterpolated<Cfg> trKrnl;
  trKrnl.tractionPlusMatrix = godunovData.tractionPlusMatrix;
  trKrnl.tractionMinusMatrix = godunovData.tractionMinusMatrix;
  trKrnl.QInterpolatedPlus = qInterpolatedPlus;
  trKrnl.QInterpolatedMinus = qInterpolatedMinus;
  trKrnl.tractionInterpolated = tractionInterpolated;
  trKrnl.execute();

  alignas(Alignment) real staticFrictionalWork[tensor::staticFrictionalWork<Cfg>::size()]{};

  dynamicRupture::kernel::accumulateStaticFrictionalWork<Cfg> feKrnl;
  feKrnl.slipInterpolated = slip;
  feKrnl.tractionInterpolated = tractionInterpolated;
  feKrnl.spaceWeights = spaceWeights;
  feKrnl.staticFrictionalWork = staticFrictionalWork;
  feKrnl.minusSurfaceArea = -0.5 * godunovData.doubledSurfaceArea;
  feKrnl.execute();

  std::array<real, multisim::NumSimulations> frictionalWorkReturn;
  std::copy(staticFrictionalWork,
            staticFrictionalWork + multisim::NumSimulations,
            std::begin(frictionalWorkReturn));
  return frictionalWorkReturn;
}

} // namespace

namespace seissol::writer {

double& EnergiesStorage::gravitationalEnergy(size_t sim) { return energies[sim][0]; }
double& EnergiesStorage::acousticEnergy(size_t sim) { return energies[sim][1]; }
double& EnergiesStorage::acousticKineticEnergy(size_t sim) { return energies[sim][2]; }
double& EnergiesStorage::elasticEnergy(size_t sim) { return energies[sim][3]; }
double& EnergiesStorage::elasticKineticEnergy(size_t sim) { return energies[sim][4]; }
double& EnergiesStorage::totalFrictionalWork(size_t sim) { return energies[sim][5]; }
double& EnergiesStorage::staticFrictionalWork(size_t sim) { return energies[sim][6]; }
double& EnergiesStorage::plasticMoment(size_t sim) { return energies[sim][7]; }
double& EnergiesStorage::seismicMoment(size_t sim) { return energies[sim][8]; }
double& EnergiesStorage::potency(size_t sim) { return energies[sim][9]; }

double& EnergiesStorage::totalMomentumX(size_t sim) { return energies[sim][10]; }
double& EnergiesStorage::totalMomentumY(size_t sim) { return energies[sim][11]; }
double& EnergiesStorage::totalMomentumZ(size_t sim) { return energies[sim][12]; }

void EnergyOutput::init(
    const GlobalData& newGlobal,
    const DynamicRupture::Storage& newDynRuptTree,
    const seissol::geometry::MeshReader& newMeshReader,
    const LTS::Storage& newStorage,
    bool newIsPlasticityEnabled,
    const std::string& outputFileNamePrefix,
    const seissol::initializer::parameters::EnergyOutputParameters& parameters) {
  if (parameters.enabled && parameters.interval > 0) {
    isEnabled = true;
  } else {
    return;
  }
  const auto rank = MPI::mpi.rank();
  logInfo() << "Initializing energy output.";

  std::size_t maxSims = 1;
  approxElements = 0;
  for (const auto& element : newMeshReader.getElements()) {
    std::visit(
        [&](auto cfg) {
          using Cfg = decltype(cfg);
          if constexpr (VolumeEnergyApproximation<model::MaterialTT<Cfg>>) {
            approxElements += 1;
          }

          maxSims = std::max(maxSims, Cfg::NumSimulations);
        },
        ConfigVariantList[element.configId]);
  }
  MPI_Allreduce(MPI_IN_PLACE,
                &approxElements,
                1,
                seissol::MPI::castToMpiType<std::size_t>(),
                MPI_SUM,
                seissol::MPI::mpi.comm());
  MPI_Allreduce(MPI_IN_PLACE,
                &maxSims,
                1,
                seissol::MPI::castToMpiType<std::size_t>(),
                MPI_MAX,
                seissol::MPI::mpi.comm());

  if (approxElements > 0) {
    logWarning() << "The volume energies printed for the given equation system are, by now, only "
                    "an isotropic approximation.";
  }

  energiesStorage.setSimcount(maxSims);
  minTimeSinceSlipRateBelowThreshold.resize(maxSims);
  minTimeSinceMomentRateBelowThreshold.resize(maxSims);
  seismicMomentPrevious.resize(maxSims);

  energyOutputInterval = parameters.interval;
  isFileOutputEnabled = rank == 0;
  isTerminalOutputEnabled = parameters.terminalOutput && (rank == 0);
  terminatorMaxTimePostRupture = parameters.terminatorMaxTimePostRupture;
  terminatorMomentRateThreshold = parameters.terminatorMomentRateThreshold;
  isCheckAbortCriteraSlipRateEnabled = std::isfinite(terminatorMaxTimePostRupture);
  isCheckAbortCriteraMomentRateEnabled = (terminatorMomentRateThreshold > 0);
  computeVolumeEnergiesEveryOutput = parameters.computeVolumeEnergiesEveryOutput;
  outputFileName = outputFileNamePrefix + "-energy.csv";

  global = &newGlobal;
  drStorage = &newDynRuptTree;
  meshReader = &newMeshReader;
  ltsStorage = &newStorage;

  isPlasticityEnabled = newIsPlasticityEnabled;

#ifdef ACL_DEVICE
  const auto maxCells = ltsStorage->getMaxClusterSize();

  if (maxCells > 0) {
    constexpr auto QSize = tensor::Q<Cfg>::size();
    timeDerivativePlusHost =
        device::DeviceInstance::getInstance().api->allocPinnedMem(maxCells * QSize * sizeof(real));
    timeDerivativeMinusHost =
        device::DeviceInstance::getInstance().api->allocPinnedMem(maxCells * QSize * sizeof(real));
    timeDerivativePlusHostMapped =
        device::DeviceInstance::getInstance().api->devicePointer(timeDerivativePlusHost);
    timeDerivativeMinusHostMapped =
        device::DeviceInstance::getInstance().api->devicePointer(timeDerivativeMinusHost);
  }
#endif

  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  setSyncInterval(parameters.interval);
}

void EnergyOutput::syncPoint(double time) {
  assert(isEnabled);
  const auto rank = MPI::mpi.rank();
  logInfo() << "Writing energy output at time" << time;
  computeEnergies();
  reduceEnergies();
  if (isCheckAbortCriteraSlipRateEnabled) {
    reduceMinTimeSinceSlipRateBelowThreshold();
  }
  if ((rank == 0) && isCheckAbortCriteraMomentRateEnabled) {
    for (size_t sim = 0; sim < energiesStorage.simcount(); sim++) {
      const double seismicMomentRate =
          (energiesStorage.seismicMoment(sim) - seismicMomentPrevious[sim]) / energyOutputInterval;
      seismicMomentPrevious[sim] = energiesStorage.seismicMoment(sim);
      if (time > 0 && seismicMomentRate < terminatorMomentRateThreshold) {
        minTimeSinceMomentRateBelowThreshold[sim] += energyOutputInterval;
      } else {
        minTimeSinceMomentRateBelowThreshold[sim] = 0.0;
      }
    }
  }
  if (isTerminalOutputEnabled) {
    printEnergies();
  }
  if (isCheckAbortCriteraSlipRateEnabled) {
    checkAbortCriterion(minTimeSinceSlipRateBelowThreshold, "All slip-rate are");
  }
  if (isCheckAbortCriteraMomentRateEnabled) {
    checkAbortCriterion(minTimeSinceMomentRateBelowThreshold, "The seismic moment rate is");
  }

  if (isFileOutputEnabled) {
    writeEnergies(time);
  }
  ++outputId;
  logInfo() << "Writing energy output at time" << time << "Done.";
}

void EnergyOutput::simulationStart(std::optional<double> checkpointTime) {
  if (isFileOutputEnabled) {
    out.open(outputFileName);
    out << std::scientific;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    writeHeader();
  }
  syncPoint(checkpointTime.value_or(0));
}

EnergyOutput::~EnergyOutput() {
#ifdef ACL_DEVICE
  if (timeDerivativePlusHost != nullptr) {
    device::DeviceInstance::getInstance().api->freePinnedMem(timeDerivativePlusHost);
  }
  if (timeDerivativeMinusHost != nullptr) {
    device::DeviceInstance::getInstance().api->freePinnedMem(timeDerivativeMinusHost);
  }
#endif
}

void EnergyOutput::computeDynamicRuptureEnergies() {
  for (size_t sim = 0; sim < energiesStorage.simcount(); sim++) {
    double& totalFrictionalWork = energiesStorage.totalFrictionalWork(sim);
    double& staticFrictionalWork = energiesStorage.staticFrictionalWork(sim);
    double& seismicMoment = energiesStorage.seismicMoment(sim);
    double& potency = energiesStorage.potency(sim);
    minTimeSinceSlipRateBelowThreshold[sim] = std::numeric_limits<double>::max();

#ifdef ACL_DEVICE
    void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();
#endif
    for (const auto& layer : drStorage->leaves()) {
      layer.wrap([&](auto cfg) {
        using Cfg = decltype(cfg);
        using real = Real<Cfg>;

        /// \todo timeDerivativePlus and timeDerivativeMinus are missing the last timestep.
        /// (We'd need to send the dofs over the network in order to fix this.)
#ifdef ACL_DEVICE
        constexpr auto QSize = tensor::Q<Cfg>::size();
        const ConditionalKey timeIntegrationKey(*KernelNames::DrTime);
        auto& table = layer.getConditionalTable<inner_keys::Dr>();
        if (table.find(timeIntegrationKey) != table.end()) {
          const auto& entry = table.at(timeIntegrationKey);
          const real** timeDerivativePlusDevice = const_cast<const real**>(
              (entry.get(inner_keys::Dr::Id::DerivativesPlus))->getDeviceDataPtr());
          const real** timeDerivativeMinusDevice = const_cast<const real**>(
              (entry.get(inner_keys::Dr::Id::DerivativesMinus))->getDeviceDataPtr());
          device::DeviceInstance::getInstance().algorithms.copyScatterToUniform(
              timeDerivativePlusDevice,
              timeDerivativePlusHostMapped,
              QSize,
              QSize,
              layer.size(),
              stream);
          device::DeviceInstance::getInstance().algorithms.copyScatterToUniform(
              timeDerivativeMinusDevice,
              timeDerivativeMinusHostMapped,
              QSize,
              QSize,
              layer.size(),
              stream);
          device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
        }
        const auto timeDerivativePlusPtr = [&](std::size_t i) {
          return timeDerivativePlusHost + QSize * i;
        };
        const auto timeDerivativeMinusPtr = [&](std::size_t i) {
          return timeDerivativeMinusHost + QSize * i;
        };
#else
        // TODO: for fused simulations, do this once and reuse
        auto* const* timeDerivativePlus = layer.var<DynamicRupture::TimeDerivativePlus>(cfg);
        auto* const* timeDerivativeMinus = layer.var<DynamicRupture::TimeDerivativeMinus>(cfg);
        const auto timeDerivativePlusPtr = [&](std::size_t i) { return timeDerivativePlus[i]; };
        const auto timeDerivativeMinusPtr = [&](std::size_t i) { return timeDerivativeMinus[i]; };
#endif
        const auto* godunovData = layer.var<DynamicRupture::GodunovData>(cfg);
        const auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
        const auto* drEnergyOutput = layer.var<DynamicRupture::DREnergyOutputVar>(cfg);
        const seissol::model::IsotropicWaveSpeeds* waveSpeedsPlus =
            layer.var<DynamicRupture::WaveSpeedsPlus>();
        const seissol::model::IsotropicWaveSpeeds* waveSpeedsMinus =
            layer.var<DynamicRupture::WaveSpeedsMinus>();
        const auto layerSize = layer.size();

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(                                                                \
        + : totalFrictionalWork, staticFrictionalWork, seismicMoment, potency)
#endif
        for (std::size_t i = 0; i < layerSize; ++i) {
          if (faceInformation[i].plusSideOnThisRank) {
            for (std::size_t j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints<Cfg>; ++j) {
              totalFrictionalWork +=
                  drEnergyOutput[i].frictionalEnergy[j * seissol::multisim::NumSimulations + sim];
            }
            staticFrictionalWork += computeStaticWork<Cfg>(timeDerivativePlusPtr(i),
                                                           timeDerivativeMinusPtr(i),
                                                           faceInformation[i],
                                                           godunovData[i],
                                                           drEnergyOutput[i].slip,
                                                           global)[sim];

            const double muPlus = waveSpeedsPlus[i].density * waveSpeedsPlus[i].sWaveVelocity *
                                  waveSpeedsPlus[i].sWaveVelocity;
            const double muMinus = waveSpeedsMinus[i].density * waveSpeedsMinus[i].sWaveVelocity *
                                   waveSpeedsMinus[i].sWaveVelocity;
            const double mu = 2.0 * muPlus * muMinus / (muPlus + muMinus);
            double potencyIncrease = 0.0;
            for (std::size_t k = 0; k < seissol::dr::misc::NumBoundaryGaussPoints<Cfg>; ++k) {
              potencyIncrease +=
                  drEnergyOutput[i].accumulatedSlip[k * seissol::multisim::NumSimulations + sim];
            }
            potencyIncrease *= 0.5 * godunovData[i].doubledSurfaceArea /
                               seissol::dr::misc::NumBoundaryGaussPoints<Cfg>;
            potency += potencyIncrease;
            seismicMoment += potencyIncrease * mu;
          }
        }
        double localMin = std::numeric_limits<double>::max();
#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(min : localMin) default(none)                                   \
    shared(layerSize, drEnergyOutput, faceInformation, sim)
#endif
        for (std::size_t i = 0; i < layerSize; ++i) {
          if (faceInformation[i].plusSideOnThisRank) {
            for (std::size_t j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints<Cfg>; ++j) {
              localMin = std::min(
                  static_cast<double>(
                      drEnergyOutput[i].timeSinceSlipRateBelowThreshold
                          [static_cast<size_t>(j * seissol::multisim::NumSimulations) + sim]),
                  localMin);
            }
          }
        }
        minTimeSinceSlipRateBelowThreshold[sim] =
            std::min(localMin, minTimeSinceSlipRateBelowThreshold[sim]);
      });
    }
  }
}

void EnergyOutput::computeVolumeEnergies() {
  for (size_t sim = 0; sim < energiesStorage.simcount(); sim++) {
    // TODO: Abstract energy calculations and implement is for anisotropic and poroelastic materials
    [[maybe_unused]] auto& totalGravitationalEnergyLocal = energiesStorage.gravitationalEnergy(sim);
    [[maybe_unused]] auto& totalAcousticEnergyLocal = energiesStorage.acousticEnergy(sim);
    [[maybe_unused]] auto& totalAcousticKineticEnergyLocal =
        energiesStorage.acousticKineticEnergy(sim);
    [[maybe_unused]] auto& totalElasticEnergyLocal = energiesStorage.elasticEnergy(sim);
    [[maybe_unused]] auto& totalElasticKineticEnergyLocal =
        energiesStorage.elasticKineticEnergy(sim);
    [[maybe_unused]] auto& totalPlasticMoment = energiesStorage.plasticMoment(sim);
    [[maybe_unused]] auto& totalMomentumXLocal = energiesStorage.totalMomentumX(sim);
    [[maybe_unused]] auto& totalMomentumYLocal = energiesStorage.totalMomentumY(sim);
    [[maybe_unused]] auto& totalMomentumZLocal = energiesStorage.totalMomentumZ(sim);

    const std::vector<Element>& elements = meshReader->getElements();
    const std::vector<Vertex>& vertices = meshReader->getVertices();

    [[maybe_unused]] const auto g = seissolInstance.getGravitationSetup().acceleration;

    // Note: Default(none) is not possible, clang requires data sharing attribute for g, gcc forbids
    // it
    for (const auto& layer : ltsStorage->leaves(Ghost)) {
      layer.wrap([&](auto cfg) {
        using Cfg = decltype(cfg);

        // NOLINTNEXTLINE
        using real = Real<Cfg>;
        using MaterialT = model::MaterialTT<Cfg>;

        const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
        const auto* cellInformationData = layer.var<LTS::CellInformation>();
        const auto* faceDisplacementsData = layer.var<LTS::FaceDisplacements>(cfg);
        const auto* materialData = layer.var<LTS::Material>();
        const auto* boundaryMappingData = layer.var<LTS::BoundaryMapping>(cfg);
        const auto* pstrainData = layer.var<LTS::PStrain>(cfg);
        const auto* dofsData = layer.var<LTS::Dofs>(cfg);
#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for schedule(static) reduction(+ : totalGravitationalEnergyLocal,             \
                                                        totalAcousticEnergyLocal,                  \
                                                        totalAcousticKineticEnergyLocal,           \
                                                        totalElasticEnergyLocal,                   \
                                                        totalElasticKineticEnergyLocal,            \
                                                        totalMomentumXLocal,                       \
                                                        totalMomentumYLocal,                       \
                                                        totalMomentumZLocal,                       \
                                                        totalPlasticMoment)                        \
    shared(elements, vertices, global)
#endif
        for (std::size_t cell = 0; cell < layer.size(); ++cell) {
          if (secondaryInformation[cell].duplicate > 0) {
            // skip duplicate cells
            continue;
          }
          const auto elementId = secondaryInformation[cell].meshId;
          const double volume = MeshTools::volume(elements[elementId], vertices);
          const auto& material = *materialData[cell].local;
          const auto& cellInformation = cellInformationData[cell];
          const auto& faceDisplacements = faceDisplacementsData[cell];

          constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
          constexpr auto NumQuadraturePointsTet = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

          double quadraturePointsTet[NumQuadraturePointsTet][3];
          double quadratureWeightsTet[NumQuadraturePointsTet];
          seissol::quadrature::TetrahedronQuadrature(
              quadraturePointsTet, quadratureWeightsTet, QuadPolyDegree);

          constexpr auto NumQuadraturePointsTri = QuadPolyDegree * QuadPolyDegree;
          double quadraturePointsTri[NumQuadraturePointsTri][2];
          double quadratureWeightsTri[NumQuadraturePointsTri];
          seissol::quadrature::TriangleQuadrature(
              quadraturePointsTri, quadratureWeightsTri, QuadPolyDegree);

          // Needed to weight the integral.
          const auto jacobiDet = 6 * volume;

          alignas(Alignment) real numericalSolutionData[tensor::dofsQP<Cfg>::size()];
          auto numericalSolution = init::dofsQP<Cfg>::view::create(numericalSolutionData);
          // Evaluate numerical solution at quad. nodes
          kernel::evalAtQP<Cfg> krnl;
          krnl.evalAtQP = global->get<Cfg>().evalAtQPMatrix;
          krnl.dofsQP = numericalSolutionData;
          krnl.Q = dofsData[cell];
          krnl.execute();

          auto numSub = multisim::simtensor(numericalSolution, sim);

          constexpr auto UIdx = MaterialT::TractionQuantities;

          for (size_t qp = 0; qp < NumQuadraturePointsTet; ++qp) {
            const auto curWeight = jacobiDet * quadratureWeightsTet[qp];
            const auto rho = material.getDensity();

            const auto u = numSub(qp, UIdx + 0);
            const auto v = numSub(qp, UIdx + 1);
            const auto w = numSub(qp, UIdx + 2);
            const double curKineticEnergy = 0.5 * rho * (u * u + v * v + w * w);
            const double curMomentumX = rho * u;
            const double curMomentumY = rho * v;
            const double curMomentumZ = rho * w;

            if (std::abs(material.getMuBar()) < 10e-14) {
              // Acoustic
              constexpr std::size_t PIdx = 0;
              const auto k = material.getLambdaBar();
              const auto p = numSub(qp, PIdx);
              const double curAcousticEnergy = (p * p) / (2 * k);
              totalAcousticEnergyLocal += curWeight * curAcousticEnergy;
              totalAcousticKineticEnergyLocal += curWeight * curKineticEnergy;
            } else {
              // Elastic
              totalElasticKineticEnergyLocal += curWeight * curKineticEnergy;
              auto getStressIndex = [](int i, int j) {
                const static auto Lookup =
                    std::array<std::array<int, 3>, 3>{{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}}};
                return Lookup[i][j];
              };
              totalMomentumXLocal += curWeight * curMomentumX;
              totalMomentumYLocal += curWeight * curMomentumY;
              totalMomentumZLocal += curWeight * curMomentumZ;

              auto getStress = [&](int i, int j) { return numSub(qp, getStressIndex(i, j)); };

              const auto lambda = material.getLambdaBar();
              const auto mu = material.getMuBar();
              const auto sumUniaxialStresses = getStress(0, 0) + getStress(1, 1) + getStress(2, 2);
              auto computeStrain = [&](int i, int j) {
                double strain = 0.0;
                const auto factor = -1.0 * (lambda) / (2.0 * mu * (3.0 * lambda + 2.0 * mu));
                if (i == j) {
                  strain += factor * sumUniaxialStresses;
                }
                strain += 1.0 / (2.0 * mu) * getStress(i, j);
                return strain;
              };
              double curElasticEnergy = 0.0;
              for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                  curElasticEnergy += getStress(i, j) * computeStrain(i, j);
                }
              }
              totalElasticEnergyLocal += curWeight * 0.5 * curElasticEnergy;
            }
          }

          const auto* boundaryMappings = boundaryMappingData[cell];
          // Compute gravitational energy
          for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
            if (cellInformation.faceTypes[face] != FaceType::FreeSurfaceGravity) {
              continue;
            }

            // Displacements are stored in face-aligned coordinate system.
            // We need to rotate it to the global coordinate system.
            const auto& boundaryMapping = boundaryMappings[face];
            auto tinv = init::Tinv<Cfg>::view::create(boundaryMapping.dataTinv);
            alignas(Alignment) real
                rotateDisplacementToFaceNormalData[init::displacementRotationMatrix<Cfg>::Size];

            auto rotateDisplacementToFaceNormal =
                init::displacementRotationMatrix<Cfg>::view::create(
                    rotateDisplacementToFaceNormalData);
            for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                rotateDisplacementToFaceNormal(i, j) = tinv(i + UIdx, j + UIdx);
              }
            }

            alignas(Alignment)
                std::array<real, tensor::rotatedFaceDisplacementAtQuadratureNodes<Cfg>::Size>
                    displQuadData{};
            const auto* curFaceDisplacementsData = faceDisplacements[face];
            seissol::kernel::rotateFaceDisplacementsAndEvaluateAtQuadratureNodes<Cfg> evalKrnl;
            evalKrnl.rotatedFaceDisplacement = curFaceDisplacementsData;
            evalKrnl.V2nTo2JacobiQuad = init::V2nTo2JacobiQuad<Cfg>::Values;
            evalKrnl.rotatedFaceDisplacementAtQuadratureNodes = displQuadData.data();
            evalKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
            evalKrnl.execute();

            // Perform quadrature
            const auto surface = MeshTools::surface(elements[elementId], face, vertices);
            const auto rho = material.getDensity();

            static_assert(NumQuadraturePointsTri ==
                          init::rotatedFaceDisplacementAtQuadratureNodes<
                              Cfg>::Shape[multisim::BasisFunctionDimension]);
            auto rotatedFaceDisplacementFused =
                init::rotatedFaceDisplacementAtQuadratureNodes<Cfg>::view::create(
                    displQuadData.data());
            auto rotatedFaceDisplacement = multisim::simtensor(rotatedFaceDisplacementFused, sim);

            for (std::size_t i = 0; i < rotatedFaceDisplacement.shape(0); ++i) {
              // See for example (Saito, Tsunami generation and propagation, 2019) section 3.2.3 for
              // derivation.
              const auto displ = rotatedFaceDisplacement(i, 0);
              const auto curEnergy = 0.5 * rho * g * displ * displ;
              const auto curWeight = 2.0 * surface * quadratureWeightsTri[i];
              totalGravitationalEnergyLocal += curWeight * curEnergy;
            }
          }

          if (isPlasticityEnabled) {
            // plastic moment
            const real* pstrainCell = pstrainData[cell];
            const double mu = material.getMuBar();
            totalPlasticMoment += mu * volume * pstrainCell[tensor::QStress<Cfg>::size() + sim];
          }
        }
      });
    }
  }
}

void EnergyOutput::computeEnergies() {
  for (auto& simEnergies : energiesStorage.energies) {
    simEnergies.fill(0.0);
  }
  if (shouldComputeVolumeEnergies()) {
    computeVolumeEnergies();
  }
  computeDynamicRuptureEnergies();
}

void EnergyOutput::reduceEnergies() {
  const auto& comm = MPI::mpi.comm();
  MPI_Allreduce(MPI_IN_PLACE,
                energiesStorage.energies.data(),
                static_cast<int>(energiesStorage.energies.size()),
                MPI_DOUBLE,
                MPI_SUM,
                comm);
}

void EnergyOutput::reduceMinTimeSinceSlipRateBelowThreshold() {
  const auto& comm = MPI::mpi.comm();
  MPI_Allreduce(MPI_IN_PLACE,
                minTimeSinceSlipRateBelowThreshold.data(),
                static_cast<int>(minTimeSinceSlipRateBelowThreshold.size()),
                MPI::castToMpiType<double>(),
                MPI_MIN,
                comm);
}

void EnergyOutput::printEnergies() {
  const auto outputPrecision =
      seissolInstance.getSeisSolParameters().output.energyParameters.terminalPrecision;

  const auto shouldPrint = [](double thresholdValue) { return std::abs(thresholdValue) > 1.e-20; };
  for (size_t sim = 0; sim < energiesStorage.simcount(); sim++) {
    const std::string fusedPrefix =
        energiesStorage.simcount() > 1 ? "[" + std::to_string(sim) + "]" : "";
    const std::string approxPrefix = approxElements > 0 ? "[approximated]" : "";
    const auto totalAcousticEnergy =
        energiesStorage.acousticKineticEnergy(sim) + energiesStorage.acousticEnergy(sim);
    const auto totalElasticEnergy =
        energiesStorage.elasticKineticEnergy(sim) + energiesStorage.elasticEnergy(sim);
    const auto ratioElasticKinematic =
        100.0 * energiesStorage.elasticKineticEnergy(sim) / totalElasticEnergy;
    const auto ratioElasticPotential =
        100.0 * energiesStorage.elasticEnergy(sim) / totalElasticEnergy;
    const auto ratioAcousticKinematic =
        100.0 * energiesStorage.acousticKineticEnergy(sim) / totalAcousticEnergy;
    const auto ratioAcousticPotential =
        100.0 * energiesStorage.acousticEnergy(sim) / totalAcousticEnergy;
    const auto totalFrictionalWork = energiesStorage.totalFrictionalWork(sim);
    const auto staticFrictionalWork = energiesStorage.staticFrictionalWork(sim);
    const auto radiatedEnergy = totalFrictionalWork - staticFrictionalWork;
    const auto ratioFrictionalStatic = 100.0 * staticFrictionalWork / totalFrictionalWork;
    const auto ratioFrictionalRadiated = 100.0 * radiatedEnergy / totalFrictionalWork;
    const auto ratioPlasticMoment =
        100.0 * energiesStorage.plasticMoment(sim) /
        (energiesStorage.plasticMoment(sim) + energiesStorage.seismicMoment(sim));
    const auto totalMomentumX = energiesStorage.totalMomentumX(sim);
    const auto totalMomentumY = energiesStorage.totalMomentumY(sim);
    const auto totalMomentumZ = energiesStorage.totalMomentumZ(sim);
    if (shouldComputeVolumeEnergies()) {
      if (shouldPrint(totalElasticEnergy)) {
        logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                  << approxPrefix.c_str()
                  << " Elastic energy (total, % kinematic, % potential): " << totalElasticEnergy
                  << " ," << ratioElasticKinematic << " ," << ratioElasticPotential;
      }
      if (shouldPrint(totalAcousticEnergy)) {
        logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                  << approxPrefix.c_str()
                  << " Acoustic energy (total, % kinematic, % potential): " << totalAcousticEnergy
                  << " ," << ratioAcousticKinematic << " ," << ratioAcousticPotential;
      }
      if (shouldPrint(energiesStorage.gravitationalEnergy(sim))) {
        logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                  << approxPrefix.c_str()
                  << " Gravitational energy:" << energiesStorage.gravitationalEnergy(sim);
      }
      if (shouldPrint(energiesStorage.plasticMoment(sim))) {
        logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                  << approxPrefix.c_str()
                  << " Plastic moment (value, equivalent Mw, % total moment):"
                  << energiesStorage.plasticMoment(sim) << " ,"
                  << 2.0 / 3.0 * std::log10(energiesStorage.plasticMoment(sim)) - 6.07 << " ,"
                  << ratioPlasticMoment;
      }
    } else {
      logInfo() << "Volume energies skipped at this step";
    }
    logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
              << " Total momentum (X, Y, Z):" << totalMomentumX << " ," << totalMomentumY << " ,"
              << totalMomentumZ;
    if (shouldPrint(totalFrictionalWork)) {
      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << " Frictional work (total, % static, % radiated): " << totalFrictionalWork << " ,"
                << ratioFrictionalStatic << " ," << ratioFrictionalRadiated;
      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << " Seismic moment (without plasticity):" << energiesStorage.seismicMoment(sim)
                << " Mw:" << 2.0 / 3.0 * std::log10(energiesStorage.seismicMoment(sim)) - 6.07;
    }
    if (!std::isfinite(totalElasticEnergy + totalAcousticEnergy)) {
      logError() << fusedPrefix << " Detected Inf/NaN in energies. Aborting.";
    }
  }
}

void EnergyOutput::checkAbortCriterion(const std::vector<double>& timeSinceThreshold,
                                       const std::string& prefixMessage) {
  size_t abortCount = 0;
  for (size_t sim = 0; sim < timeSinceThreshold.size(); sim++) {
    if ((timeSinceThreshold[sim] > 0) and
        (timeSinceThreshold[sim] < std::numeric_limits<double>::max())) {
      if (static_cast<double>(timeSinceThreshold[sim]) < terminatorMaxTimePostRupture) {
        logInfo() << prefixMessage.c_str() << "below threshold since" << timeSinceThreshold[sim]
                  << "in simulation: " << sim
                  << "s (lower than the abort criteria: " << terminatorMaxTimePostRupture << "s)";
      } else {
        logInfo() << prefixMessage.c_str() << "below threshold since" << timeSinceThreshold[sim]
                  << "in simulation: " << sim
                  << "s (greater than the abort criteria: " << terminatorMaxTimePostRupture << "s)";
        ++abortCount;
      }
    }
  }

  bool abort = abortCount == timeSinceThreshold.size();
  const auto& comm = MPI::mpi.comm();
  MPI_Bcast(reinterpret_cast<void*>(&abort), 1, MPI_CXX_BOOL, 0, comm);
  if (abort) {
    seissolInstance.simulator().abort();
  }
}

void EnergyOutput::writeHeader() {
  out << "time,variable,simulation_index,measurement" << std::endl;
}

void EnergyOutput::writeEnergies(double time) {
  for (size_t sim = 0; sim < energiesStorage.simcount(); sim++) {
    const std::string fusedSuffix = std::to_string(sim);
    if (shouldComputeVolumeEnergies()) {
      out << time << ",gravitational_energy," << fusedSuffix << ","
          << energiesStorage.gravitationalEnergy(sim) << "\n"
          << time << ",acoustic_energy," << fusedSuffix << ","
          << energiesStorage.acousticEnergy(sim) << "\n"
          << time << ",acoustic_kinetic_energy," << fusedSuffix << ","
          << energiesStorage.acousticKineticEnergy(sim) << "\n"
          << time << ",elastic_energy," << fusedSuffix << "," << energiesStorage.elasticEnergy(sim)
          << "\n"
          << time << ",elastic_kinetic_energy," << fusedSuffix << ","
          << energiesStorage.elasticKineticEnergy(sim) << "\n"
          << time << ",plastic_moment," << fusedSuffix << "," << energiesStorage.plasticMoment(sim)
          << "\n"
          << time << ",momentumX," << fusedSuffix << "," << energiesStorage.totalMomentumX(sim)
          << "\n"
          << time << ",momentumY," << fusedSuffix << "," << energiesStorage.totalMomentumY(sim)
          << "\n"
          << time << ",momentumZ," << fusedSuffix << "," << energiesStorage.totalMomentumZ(sim)
          << "\n";
    }
    out << time << ",total_frictional_work," << fusedSuffix << ","
        << energiesStorage.totalFrictionalWork(sim) << "\n"
        << time << ",static_frictional_work," << fusedSuffix << ","
        << energiesStorage.staticFrictionalWork(sim) << "\n"
        << time << ",seismic_moment," << fusedSuffix << "," << energiesStorage.seismicMoment(sim)
        << "\n"
        << time << ",potency," << fusedSuffix << "," << energiesStorage.potency(sim) << "\n"
        << time << ",plastic_moment," << fusedSuffix << "," << energiesStorage.plasticMoment(sim)
        << std::endl;
  }
}

bool EnergyOutput::shouldComputeVolumeEnergies() const {
  return outputId % computeVolumeEnergiesEveryOutput == 0;
}

} // namespace seissol::writer
