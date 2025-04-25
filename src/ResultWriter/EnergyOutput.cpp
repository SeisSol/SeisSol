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
#include <Common/Constants.h>
#include <Equations/Datastructures.h>
#include <Geometry/MeshDefinition.h>
#include <Geometry/MeshTools.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Parameters/OutputParameters.h>
#include <Initializer/PreProcessorMacros.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Lut.h>
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
#include <init.h>
#include <iomanip>
#include <ios>
#include <kernel.h>
#include <limits>
#include <mpi.h>
#include <optional>
#include <ostream>
#include <string>
#include <tensor.h>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#endif

namespace {
constexpr bool VolumeEnergyApproximation =
    model::MaterialT::Type != model::MaterialType::Elastic &&
    model::MaterialT::Type != model::MaterialType::Viscoelastic &&
    model::MaterialT::Type != model::MaterialType::Acoustic;
} // namespace

namespace seissol::writer {

double& EnergiesStorage::gravitationalEnergy(size_t sim) {
  return energies[0 + sim * NumberOfEnergies];
}
double& EnergiesStorage::acousticEnergy(size_t sim) { return energies[1 + sim * NumberOfEnergies]; }
double& EnergiesStorage::acousticKineticEnergy(size_t sim) {
  return energies[2 + sim * NumberOfEnergies];
}
double& EnergiesStorage::elasticEnergy(size_t sim) { return energies[3 + sim * NumberOfEnergies]; }
double& EnergiesStorage::elasticKineticEnergy(size_t sim) {
  return energies[4 + sim * NumberOfEnergies];
}
double& EnergiesStorage::totalFrictionalWork(size_t sim) {
  return energies[5 + sim * NumberOfEnergies];
}
double& EnergiesStorage::staticFrictionalWork(size_t sim) {
  return energies[6 + sim * NumberOfEnergies];
}
double& EnergiesStorage::plasticMoment(size_t sim) { return energies[7 + sim * NumberOfEnergies]; }
double& EnergiesStorage::seismicMoment(size_t sim) { return energies[8 + sim * NumberOfEnergies]; }
double& EnergiesStorage::potency(size_t sim) { return energies[9 + sim * NumberOfEnergies]; }

double& EnergiesStorage::totalMomentumX(size_t sim) {
  return energies[10 + sim * NumberOfEnergies];
}
double& EnergiesStorage::totalMomentumY(size_t sim) {
  return energies[11 + sim * NumberOfEnergies];
}
double& EnergiesStorage::totalMomentumZ(size_t sim) {
  return energies[12 + sim * NumberOfEnergies];
}

void EnergyOutput::init(
    GlobalData* newGlobal,
    seissol::initializer::DynamicRupture* newDynRup,
    seissol::initializer::LTSTree* newDynRuptTree,
    seissol::geometry::MeshReader* newMeshReader,
    seissol::initializer::LTSTree* newLtsTree,
    seissol::initializer::LTS* newLts,
    seissol::initializer::Lut* newLtsLut,
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

  if constexpr (VolumeEnergyApproximation) {
    logWarning() << "The volume energies printed for the given equation system"
                 << model::MaterialT::Text << "are, by now, only an isotropic approximation.";
  }

  energyOutputInterval = parameters.interval;
  isFileOutputEnabled = rank == 0;
  isTerminalOutputEnabled = parameters.terminalOutput && (rank == 0);
  terminatorMaxTimePostRupture = parameters.terminatorMaxTimePostRupture;
  terminatorMomentRateThreshold = parameters.terminatorMomentRateThreshold;
  isCheckAbortCriteraSlipRateEnabled = std::isfinite(terminatorMaxTimePostRupture);
  isCheckAbortCriteraMomentRateEnabled = (terminatorMomentRateThreshold > 0);
  computeVolumeEnergiesEveryOutput = parameters.computeVolumeEnergiesEveryOutput;
  outputFileName = outputFileNamePrefix + "-energy.csv";

  global = newGlobal;
  dynRup = newDynRup;
  dynRupTree = newDynRuptTree;
  meshReader = newMeshReader;
  ltsTree = newLtsTree;
  lts = newLts;
  ltsLut = newLtsLut;

  isPlasticityEnabled = newIsPlasticityEnabled;

#ifdef ACL_DEVICE
  const auto maxCells = ltsTree->getMaxClusterSize();

  if (maxCells > 0) {
    constexpr auto QSize = tensor::Q::size();
    timeDerivativePlusHost = reinterpret_cast<real*>(
        device::DeviceInstance::getInstance().api->allocPinnedMem(maxCells * QSize * sizeof(real)));
    timeDerivativeMinusHost = reinterpret_cast<real*>(
        device::DeviceInstance::getInstance().api->allocPinnedMem(maxCells * QSize * sizeof(real)));
    timeDerivativePlusHostMapped = reinterpret_cast<real*>(
        device::DeviceInstance::getInstance().api->devicePointer(timeDerivativePlusHost));
    timeDerivativeMinusHostMapped = reinterpret_cast<real*>(
        device::DeviceInstance::getInstance().api->devicePointer(timeDerivativeMinusHost));
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
    for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
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
    out << std::setprecision(std::numeric_limits<real>::max_digits10);
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

std::array<real, multisim::NumSimulations>
    EnergyOutput::computeStaticWork(const real* degreesOfFreedomPlus,
                                    const real* degreesOfFreedomMinus,
                                    const DRFaceInformation& faceInfo,
                                    const DRGodunovData& godunovData,
                                    const real slip[seissol::tensor::slipInterpolated::size()]) {
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
  krnl.TinvT = godunovData.TinvT;
  krnl._prefetch.QInterpolated = qInterpolatedPlus;
  krnl.execute(faceInfo.plusSide, 0);

  krnl.QInterpolated = qInterpolatedMinus;
  krnl.Q = qMinus;
  krnl.TinvT = godunovData.TinvT;
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

  std::array<real, multisim::NumSimulations> frictionalWorkReturn;
  std::copy(staticFrictionalWork,
            staticFrictionalWork + multisim::NumSimulations,
            std::begin(frictionalWorkReturn));
  return frictionalWorkReturn;
}

void EnergyOutput::computeDynamicRuptureEnergies() {
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    double& totalFrictionalWork = energiesStorage.totalFrictionalWork(sim);
    double& staticFrictionalWork = energiesStorage.staticFrictionalWork(sim);
    double& seismicMoment = energiesStorage.seismicMoment(sim);
    double& potency = energiesStorage.potency(sim);
    minTimeSinceSlipRateBelowThreshold[sim] = std::numeric_limits<real>::max();

#ifdef ACL_DEVICE
    void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();
#endif
    for (auto& layer : dynRupTree->leaves()) {
      /// \todo timeDerivativePlus and timeDerivativeMinus are missing the last timestep.
      /// (We'd need to send the dofs over the network in order to fix this.)
#ifdef ACL_DEVICE
      constexpr auto QSize = tensor::Q::size();
      const ConditionalKey timeIntegrationKey(*KernelNames::DrTime);
      auto& table = layer.getConditionalTable<inner_keys::Dr>();
      if (table.find(timeIntegrationKey) != table.end()) {
        auto& entry = table[timeIntegrationKey];
        real** timeDerivativePlusDevice =
            (entry.get(inner_keys::Dr::Id::DerivativesPlus))->getDeviceDataPtr();
        real** timeDerivativeMinusDevice =
            (entry.get(inner_keys::Dr::Id::DerivativesMinus))->getDeviceDataPtr();
        device::DeviceInstance::getInstance().algorithms.copyScatterToUniform(
            timeDerivativePlusDevice,
            timeDerivativePlusHostMapped,
            QSize,
            QSize,
            layer.getNumberOfCells(),
            stream);
        device::DeviceInstance::getInstance().algorithms.copyScatterToUniform(
            timeDerivativeMinusDevice,
            timeDerivativeMinusHostMapped,
            QSize,
            QSize,
            layer.getNumberOfCells(),
            stream);
        device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
      }
      const auto timeDerivativePlusPtr = [&](unsigned i) {
        return timeDerivativePlusHost + QSize * i;
      };
      const auto timeDerivativeMinusPtr = [&](unsigned i) {
        return timeDerivativeMinusHost + QSize * i;
      };
#else
      real** timeDerivativePlus = layer.var(dynRup->timeDerivativePlus);
      real** timeDerivativeMinus = layer.var(dynRup->timeDerivativeMinus);
      const auto timeDerivativePlusPtr = [&](unsigned i) { return timeDerivativePlus[i]; };
      const auto timeDerivativeMinusPtr = [&](unsigned i) { return timeDerivativeMinus[i]; };
#endif
      DRGodunovData* godunovData = layer.var(dynRup->godunovData);
      DRFaceInformation* faceInformation = layer.var(dynRup->faceInformation);
      DREnergyOutput* drEnergyOutput = layer.var(dynRup->drEnergyOutput);
      seissol::model::IsotropicWaveSpeeds* waveSpeedsPlus = layer.var(dynRup->waveSpeedsPlus);
      seissol::model::IsotropicWaveSpeeds* waveSpeedsMinus = layer.var(dynRup->waveSpeedsMinus);
      const auto layerSize = layer.getNumberOfCells();

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(                                                                \
        + : totalFrictionalWork, staticFrictionalWork, seismicMoment, potency) default(none)       \
    shared(layerSize,                                                                              \
               drEnergyOutput,                                                                     \
               faceInformation,                                                                    \
               timeDerivativeMinusPtr,                                                             \
               timeDerivativePlusPtr,                                                              \
               godunovData,                                                                        \
               waveSpeedsPlus,                                                                     \
               waveSpeedsMinus,                                                                    \
               sim)
#endif
      for (unsigned i = 0; i < layerSize; ++i) {
        if (faceInformation[i].plusSideOnThisRank) {
          for (unsigned j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints; ++j) {
            totalFrictionalWork +=
                drEnergyOutput[i].frictionalEnergy[j * seissol::multisim::NumSimulations + sim];
          }
          staticFrictionalWork += computeStaticWork(timeDerivativePlusPtr(i),
                                                    timeDerivativeMinusPtr(i),
                                                    faceInformation[i],
                                                    godunovData[i],
                                                    drEnergyOutput[i].slip)[sim];

          const real muPlus = waveSpeedsPlus[i].density * waveSpeedsPlus[i].sWaveVelocity *
                              waveSpeedsPlus[i].sWaveVelocity;
          const real muMinus = waveSpeedsMinus[i].density * waveSpeedsMinus[i].sWaveVelocity *
                               waveSpeedsMinus[i].sWaveVelocity;
          const real mu = 2.0 * muPlus * muMinus / (muPlus + muMinus);
          real potencyIncrease = 0.0;
          for (unsigned k = 0; k < seissol::dr::misc::NumBoundaryGaussPoints; ++k) {
            potencyIncrease +=
                drEnergyOutput[i].accumulatedSlip[k * seissol::multisim::NumSimulations + sim];
          }
          potencyIncrease *=
              0.5 * godunovData[i].doubledSurfaceArea / seissol::dr::misc::NumBoundaryGaussPoints;
          potency += potencyIncrease;
          seismicMoment += potencyIncrease * mu;
        }
      }
      real localMin = std::numeric_limits<real>::max();
#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(min : localMin) default(none)                                   \
    shared(layerSize, drEnergyOutput, faceInformation, sim)
#endif
      for (unsigned i = 0; i < layerSize; ++i) {
        if (faceInformation[i].plusSideOnThisRank) {
          for (unsigned j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints; ++j) {
            localMin = std::min(
                drEnergyOutput[i]
                    .timeSinceSlipRateBelowThreshold[j * seissol::multisim::NumSimulations + sim],
                localMin);
          }
        }
      }
      minTimeSinceSlipRateBelowThreshold[sim] =
          std::min(localMin, minTimeSinceSlipRateBelowThreshold[sim]);
    }
  }
}

void EnergyOutput::computeVolumeEnergies() {
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
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
    shared(elements, vertices, lts, ltsLut, global)
#endif
    for (std::size_t elementId = 0; elementId < elements.size(); ++elementId) {
      const real volume = MeshTools::volume(elements[elementId], vertices);
      const CellMaterialData& material = ltsLut->lookup(lts->material, elementId);
      auto& cellInformation = ltsLut->lookup(lts->cellInformation, elementId);
      auto& faceDisplacements = ltsLut->lookup(lts->faceDisplacements, elementId);

      constexpr auto QuadPolyDegree = ConvergenceOrder + 1;
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

      alignas(Alignment) real numericalSolutionData[tensor::dofsQP::size()];
      auto numericalSolution = init::dofsQP::view::create(numericalSolutionData);
      // Evaluate numerical solution at quad. nodes
      kernel::evalAtQP krnl;
      krnl.evalAtQP = global->evalAtQPMatrix;
      krnl.dofsQP = numericalSolutionData;
      krnl.Q = ltsLut->lookup(lts->dofs, elementId);
      krnl.execute();

      auto numSub = multisim::simtensor(numericalSolution, sim);

      // TODO: move to the material class (maybe done by #1297 + MaterialT::NumTractionQuantities)
      constexpr int UIdx = model::MaterialT::Type == model::MaterialType::Acoustic ? 1 : 6;

      for (size_t qp = 0; qp < NumQuadraturePointsTet; ++qp) {
        const auto curWeight = jacobiDet * quadratureWeightsTet[qp];
        const auto rho = material.local.rho;

        const auto u = numSub(qp, UIdx + 0);
        const auto v = numSub(qp, UIdx + 1);
        const auto w = numSub(qp, UIdx + 2);
        const double curKineticEnergy = 0.5 * rho * (u * u + v * v + w * w);
        const double curMomentumX = rho * u;
        const double curMomentumY = rho * v;
        const double curMomentumZ = rho * w;

        if (std::abs(material.local.getMuBar()) < 10e-14) {
          // Acoustic
          constexpr int PIdx = 0;
          const auto k = material.local.getLambdaBar();
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

          const auto lambda = material.local.getLambdaBar();
          const auto mu = material.local.getMuBar();
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

      auto* boundaryMappings = ltsLut->lookup(lts->boundaryMapping, elementId);
      // Compute gravitational energy
      for (int face = 0; face < 4; ++face) {
        if (cellInformation.faceTypes[face] != FaceType::FreeSurfaceGravity) {
          continue;
        }

        // Displacements are stored in face-aligned coordinate system.
        // We need to rotate it to the global coordinate system.
        auto& boundaryMapping = boundaryMappings[face];
        auto tinv = init::Tinv::view::create(boundaryMapping.TinvData);
        alignas(Alignment)
            real rotateDisplacementToFaceNormalData[init::displacementRotationMatrix::Size];

        auto rotateDisplacementToFaceNormal =
            init::displacementRotationMatrix::view::create(rotateDisplacementToFaceNormalData);
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            rotateDisplacementToFaceNormal(i, j) = tinv(i + UIdx, j + UIdx);
          }
        }

        alignas(Alignment) std::array<real, tensor::rotatedFaceDisplacementAtQuadratureNodes::Size>
            displQuadData{};
        const auto* curFaceDisplacementsData = faceDisplacements[face];
        seissol::kernel::rotateFaceDisplacementsAndEvaluateAtQuadratureNodes evalKrnl;
        evalKrnl.rotatedFaceDisplacement = curFaceDisplacementsData;
        evalKrnl.V2nTo2JacobiQuad = init::V2nTo2JacobiQuad::Values;
        evalKrnl.rotatedFaceDisplacementAtQuadratureNodes = displQuadData.data();
        evalKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
        evalKrnl.execute();

        // Perform quadrature
        const auto surface = MeshTools::surface(elements[elementId], face, vertices);
        const auto rho = material.local.rho;

        static_assert(NumQuadraturePointsTri == init::rotatedFaceDisplacementAtQuadratureNodes::
                                                    Shape[multisim::BasisFunctionDimension]);
        auto rotatedFaceDisplacementFused =
            init::rotatedFaceDisplacementAtQuadratureNodes::view::create(displQuadData.data());
        auto rotatedFaceDisplacement = multisim::simtensor(rotatedFaceDisplacementFused, sim);

        for (unsigned i = 0; i < rotatedFaceDisplacement.shape(0); ++i) {
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
        real* pstrainCell = ltsLut->lookup(lts->pstrain, elementId);
        const real mu = material.local.getMuBar();
        totalPlasticMoment += mu * volume * pstrainCell[tensor::QStress::size() + sim];
      }
    }
  }
}

void EnergyOutput::computeEnergies() {
  energiesStorage.energies.fill(0.0);
  if (shouldComputeVolumeEnergies()) {
    computeVolumeEnergies();
  }
  computeDynamicRuptureEnergies();
}

void EnergyOutput::reduceEnergies() {
#ifdef USE_MPI
  const auto rank = MPI::mpi.rank();
  const auto& comm = MPI::mpi.comm();
  MPI_Allreduce(MPI_IN_PLACE,
                energiesStorage.energies.data(),
                static_cast<int>(energiesStorage.energies.size()),
                MPI_DOUBLE,
                MPI_SUM,
                comm);
#endif
}

void EnergyOutput::reduceMinTimeSinceSlipRateBelowThreshold() {
#ifdef USE_MPI
  const auto rank = MPI::mpi.rank();
  const auto& comm = MPI::mpi.comm();
  MPI_Allreduce(MPI_IN_PLACE,
                minTimeSinceSlipRateBelowThreshold.data(),
                static_cast<int>(minTimeSinceSlipRateBelowThreshold.size()),
                MPI::castToMpiType<real>(),
                MPI_MIN,
                comm);
#endif
}

void EnergyOutput::printEnergies() {
  const auto rank = MPI::mpi.rank();
  const auto outputPrecision =
      seissolInstance.getSeisSolParameters().output.energyParameters.terminalPrecision;

  const auto shouldPrint = [](double thresholdValue) { return std::abs(thresholdValue) > 1.e-20; };
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    const std::string fusedPrefix =
        multisim::MultisimEnabled ? "[" + std::to_string(sim) + "]" : "";
    const std::string approxPrefix = VolumeEnergyApproximation ? "[approximated]" : "";
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
    if (shouldPrint(totalAcousticEnergy)) {
      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << " Acoustic energy (total, % kinematic, % potential): " << totalAcousticEnergy
                << " ," << ratioAcousticKinematic << " ," << ratioAcousticPotential;
    }
    if (shouldPrint(energiesStorage.gravitationalEnergy(sim))) {
      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << " Gravitational energy:" << energiesStorage.gravitationalEnergy(sim);
    }
    if (shouldPrint(energiesStorage.plasticMoment(sim))) {
      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << " Plastic moment (value, equivalent Mw, % total moment):"
                << energiesStorage.plasticMoment(sim) << " ,"
                << 2.0 / 3.0 * std::log10(energiesStorage.plasticMoment(sim)) - 6.07 << " ,"
                << ratioPlasticMoment;
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

void EnergyOutput::checkAbortCriterion(
    const std::array<real, multisim::NumSimulations>& timeSinceThreshold,
    const std::string& prefixMessage) {
  const auto rank = MPI::mpi.rank();
  bool abort = true;
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    if ((timeSinceThreshold[sim] > 0) and
        (timeSinceThreshold[sim] < std::numeric_limits<real>::max())) {
      if (static_cast<double>(timeSinceThreshold[sim]) < terminatorMaxTimePostRupture) {
        logInfo() << prefixMessage.c_str() << "below threshold since" << timeSinceThreshold[sim]
                  << "in simulation: " << sim
                  << "s (lower than the abort criteria: " << terminatorMaxTimePostRupture << "s)";
        abort = false;
      } else {
        logInfo() << prefixMessage.c_str() << "below threshold since" << timeSinceThreshold[sim]
                  << "in simulation: " << sim
                  << "s (greater than the abort criteria: " << terminatorMaxTimePostRupture << "s)";
      }
    }
  }
#ifdef USE_MPI
  const auto& comm = MPI::mpi.comm();
  MPI_Bcast(reinterpret_cast<void*>(&abort), 1, MPI_CXX_BOOL, 0, comm);
#endif
  if (abort) {
    seissolInstance.simulator().abort();
  }
}

void EnergyOutput::writeHeader() {
  out << "time,variable,simulation_index,measurement" << std::endl;
}

void EnergyOutput::writeEnergies(double time) {
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
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
