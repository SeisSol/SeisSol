#include "EnergyOutput.h"

#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Kernels/DynamicRupture.h"
#include "Numerical_aux/Quadrature.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include <array>

namespace seissol { namespace writer {

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

double& EnergiesStorage::totalMomentumX() { return energies[10]; }
double& EnergiesStorage::totalMomentumY() { return energies[11]; }
double& EnergiesStorage::totalMomentumZ() { return energies[12]; }

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
  logInfo(rank) << "Initializing energy output.";

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

  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  setSyncInterval(parameters.interval);
}

void EnergyOutput::syncPoint(double time) {
  assert(isEnabled);
  const auto rank = MPI::mpi.rank();
  logInfo(rank) << "Writing energy output at time" << time;
  computeEnergies();
  reduceEnergies();
  if (isCheckAbortCriteraSlipRateEnabled) {
    reduceMinTimeSinceSlipRateBelowThreshold();
  }
  if ((rank == 0) && isCheckAbortCriteraMomentRateEnabled) {
    double seismicMomentRate =
        (energiesStorage.seismicMoment() - seismicMomentPrevious) / energyOutputInterval;
    seismicMomentPrevious = energiesStorage.seismicMoment();
    if (time > 0 && seismicMomentRate < terminatorMomentRateThreshold) {
      minTimeSinceMomentRateBelowThreshold += energyOutputInterval;
    } else {
      minTimeSinceMomentRateBelowThreshold = 0.0;
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
  logInfo(rank) << "Writing energy output at time" << time << "Done.";
}

void EnergyOutput::simulationStart() {
  if (isFileOutputEnabled) {
    out.open(outputFileName);
    writeHeader();
  }
  syncPoint(0.0);
}

std::array<real, multipleSimulations::numberOfSimulations>
    EnergyOutput::computeStaticWork(const real* degreesOfFreedomPlus,
                                    const real* degreesOfFreedomMinus,
                                    const DRFaceInformation& faceInfo,
                                    const DRGodunovData& godunovData,
                                    const real slip[seissol::tensor::slipInterpolated::size()]) {
  real points[NUMBER_OF_SPACE_QUADRATURE_POINTS][2];
  alignas(ALIGNMENT) real spaceWeights[NUMBER_OF_SPACE_QUADRATURE_POINTS];
  seissol::quadrature::TriangleQuadrature(points, spaceWeights, CONVERGENCE_ORDER + 1);

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints krnl;
  krnl.V3mTo2n = global->faceToNodalMatrices;

  alignas(PAGESIZE_STACK) real QInterpolatedPlus[tensor::QInterpolatedPlus::size()];
  alignas(PAGESIZE_STACK) real QInterpolatedMinus[tensor::QInterpolatedMinus::size()];
  alignas(ALIGNMENT) real tractionInterpolated[tensor::tractionInterpolated::size()];
  alignas(ALIGNMENT) real QPlus[tensor::Q::size()];
  alignas(ALIGNMENT) real QMinus[tensor::Q::size()];

  // needed to counter potential mis-alignment
  std::memcpy(QPlus, degreesOfFreedomPlus, sizeof(QPlus));
  std::memcpy(QMinus, degreesOfFreedomMinus, sizeof(QMinus));

  krnl.QInterpolated = QInterpolatedPlus;
  krnl.Q = QPlus;
  krnl.TinvT = godunovData.TinvT;
  krnl._prefetch.QInterpolated = QInterpolatedPlus;
  krnl.execute(faceInfo.plusSide, 0);

  krnl.QInterpolated = QInterpolatedMinus;
  krnl.Q = QMinus;
  krnl.TinvT = godunovData.TinvT;
  krnl._prefetch.QInterpolated = QInterpolatedMinus;
  krnl.execute(faceInfo.minusSide, faceInfo.faceRelation);

  dynamicRupture::kernel::computeTractionInterpolated trKrnl;
  trKrnl.tractionPlusMatrix = godunovData.tractionPlusMatrix;
  trKrnl.tractionMinusMatrix = godunovData.tractionMinusMatrix;
  trKrnl.QInterpolatedPlus = QInterpolatedPlus;
  trKrnl.QInterpolatedMinus = QInterpolatedMinus;
  trKrnl.tractionInterpolated = tractionInterpolated;
  trKrnl.execute();

  alignas(ALIGNMENT) real staticFrictionalWork[tensor::staticFrictionalWork::size()]{};

  dynamicRupture::kernel::accumulateStaticFrictionalWork feKrnl;
  feKrnl.slipInterpolated = slip;
  feKrnl.tractionInterpolated = tractionInterpolated;
  feKrnl.spaceWeights = spaceWeights;
  feKrnl.staticFrictionalWork = staticFrictionalWork;
  feKrnl.minusSurfaceArea = -0.5 * godunovData.doubledSurfaceArea;
  feKrnl.execute();

  std::array<real, multipleSimulations::numberOfSimulations> frictionalWorkReturn;
  std::copy(staticFrictionalWork,
            staticFrictionalWork + multipleSimulations::numberOfSimulations,
            std::begin(frictionalWorkReturn));
  return frictionalWorkReturn;
}

void EnergyOutput::computeDynamicRuptureEnergies() {
  for (size_t sim = 0; sim < multipleSimulations::numberOfSimulations; sim++) {
    double& totalFrictionalWork = energiesStorage.totalFrictionalWork(sim);
    double& staticFrictionalWork = energiesStorage.staticFrictionalWork(sim);
    double& seismicMoment = energiesStorage.seismicMoment(sim);
    double& potency = energiesStorage.potency(sim);
    minTimeSinceSlipRateBelowThreshold = std::numeric_limits<real>::max();

#ifdef ACL_DEVICE
    unsigned maxCells = 0;
    for (auto it = dynRupTree->beginLeaf(); it != dynRupTree->endLeaf(); ++it) {
      maxCells = std::max(it->getNumberOfCells(), maxCells);
    }

    void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();

    constexpr auto qSize = tensor::Q::size();
    real* timeDerivativePlusHost = reinterpret_cast<real*>(
        device::DeviceInstance::getInstance().api->allocPinnedMem(maxCells * qSize * sizeof(real)));
    real* timeDerivativeMinusHost = reinterpret_cast<real*>(
        device::DeviceInstance::getInstance().api->allocPinnedMem(maxCells * qSize * sizeof(real)));
    const auto timeDerivativePlusPtr = [&](unsigned i) {
      return timeDerivativePlusHost + qSize * i;
    };
    const auto timeDerivativeMinusPtr = [&](unsigned i) {
      return timeDerivativeMinusHost + qSize * i;
    };
#else
    real** timeDerivativePlus = it->var(dynRup->timeDerivativePlus);
    real** timeDerivativeMinus = it->var(dynRup->timeDerivativeMinus);
    const auto timeDerivativePlusPtr = [&](unsigned i) { return timeDerivativePlus[i]; };
    const auto timeDerivativeMinusPtr = [&](unsigned i) { return timeDerivativeMinus[i]; };
#endif
    for (auto it = dynRupTree->beginLeaf(); it != dynRupTree->endLeaf(); ++it) {
      /// \todo timeDerivativePlus and timeDerivativeMinus are missing the last timestep.
      /// (We'd need to send the dofs over the network in order to fix this.)
#ifdef ACL_DEVICE
      ConditionalKey timeIntegrationKey(*KernelNames::DrTime);
      auto& table = it->getConditionalTable<inner_keys::Dr>();
      if (table.find(timeIntegrationKey) != table.end()) {
        auto& entry = table[timeIntegrationKey];
        real** timeDerivativePlusDevice =
            (entry.get(inner_keys::Dr::Id::DerivativesPlus))->getDeviceDataPtr();
        real** timeDerivativeMinusDevice =
            (entry.get(inner_keys::Dr::Id::DerivativesMinus))->getDeviceDataPtr();
        device::DeviceInstance::getInstance().algorithms.copyScatterToUniform(
            timeDerivativePlusDevice,
            timeDerivativePlusHost,
            qSize,
            qSize,
            it->getNumberOfCells(),
            stream);
        device::DeviceInstance::getInstance().algorithms.copyScatterToUniform(
            timeDerivativeMinusDevice,
            timeDerivativeMinusHost,
            qSize,
            qSize,
            it->getNumberOfCells(),
            stream);
        device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
      }
      auto const timeDerivativePlusPtr = [&](unsigned i) {
        return timeDerivativePlusHost + qSize * i;
      };
      auto const timeDerivativeMinusPtr = [&](unsigned i) {
        return timeDerivativeMinusHost + qSize * i;
      };
#else
      real** timeDerivativePlus = it->var(dynRup->timeDerivativePlus);
      real** timeDerivativeMinus = it->var(dynRup->timeDerivativeMinus);
      auto const timeDerivativePlusPtr = [&](unsigned i) { return timeDerivativePlus[i]; };
      auto const timeDerivativeMinusPtr = [&](unsigned i) { return timeDerivativeMinus[i]; };
#endif
      DRGodunovData* godunovData = it->var(dynRup->godunovData);
      DRFaceInformation* faceInformation = it->var(dynRup->faceInformation);
      DREnergyOutput* drEnergyOutput = it->var(dynRup->drEnergyOutput);
      seissol::model::IsotropicWaveSpeeds* waveSpeedsPlus = it->var(dynRup->waveSpeedsPlus);
      seissol::model::IsotropicWaveSpeeds* waveSpeedsMinus = it->var(dynRup->waveSpeedsMinus);

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(                                                                \
        + : totalFrictionalWork, staticFrictionalWork, seismicMoment, potency) default(none)       \
    shared(it,                                                                                     \
               drEnergyOutput,                                                                     \
               faceInformation,                                                                    \
               timeDerivativeMinusPtr,                                                             \
               timeDerivativePlusPtr,                                                              \
               godunovData,                                                                        \
               waveSpeedsPlus,                                                                     \
               waveSpeedsMinus,                                                                    \
               sim)
#endif
      for (unsigned i = 0; i < it->getNumberOfCells(); ++i) {
        if (faceInformation[i].plusSideOnThisRank) {
          for (unsigned j = 0; j < seissol::dr::misc::numberOfBoundaryGaussPoints; ++j) {
            totalFrictionalWork += drEnergyOutput[i].frictionalEnergy[j];
          }
          // TODO: move this out of the sim loop
          staticFrictionalWork += computeStaticWork(timeDerivativePlusPtr(i),
                                                    timeDerivativeMinusPtr(i),
                                                    faceInformation[i],
                                                    godunovData[i],
                                                    drEnergyOutput[i].slip)[sim];

          real muPlus = waveSpeedsPlus[i].density * waveSpeedsPlus[i].sWaveVelocity *
                        waveSpeedsPlus[i].sWaveVelocity;
          real muMinus = waveSpeedsMinus[i].density * waveSpeedsMinus[i].sWaveVelocity *
                         waveSpeedsMinus[i].sWaveVelocity;
          real mu = 2.0 * muPlus * muMinus / (muPlus + muMinus);
          real potencyIncrease = 0.0;
          for (unsigned k = 0; k < seissol::dr::misc::numberOfBoundaryGaussPoints; ++k) {
            potencyIncrease += drEnergyOutput[i].accumulatedSlip[k];
          }
          potencyIncrease *= 0.5 * godunovData[i].doubledSurfaceArea /
                             seissol::dr::misc::numberOfBoundaryGaussPoints;
          potency += potencyIncrease;
          seismicMoment += potencyIncrease * mu;
        }
      }
      real localMin = std::numeric_limits<real>::max();

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(min : localMin) default(none)                                   \
    shared(it, drEnergyOutput, faceInformation, minTimeSinceSlipRateBelowThreshold)
#endif
      for (unsigned i = 0; i < it->getNumberOfCells(); ++i) {
        if (faceInformation[i].plusSideOnThisRank) {
          for (unsigned j = 0; j < seissol::dr::misc::numberOfBoundaryGaussPoints; ++j) {
            if (drEnergyOutput[i].timeSinceSlipRateBelowThreshold[j] < localMin) {
              localMin = drEnergyOutput[i].timeSinceSlipRateBelowThreshold[j];
            }
          }
        }
      }
      minTimeSinceSlipRateBelowThreshold = std::min(localMin, minTimeSinceSlipRateBelowThreshold);
    }
  }
#ifdef ACL_DEVICE
  device::DeviceInstance::getInstance().api->freePinnedMem(timeDerivativePlusHost);
  device::DeviceInstance::getInstance().api->freePinnedMem(timeDerivativeMinusHost);
#endif
}

void EnergyOutput::computeVolumeEnergies() {
  for (size_t sim = 0; sim < multipleSimulations::numberOfSimulations; sim++) {
    auto& totalGravitationalEnergyLocal = energiesStorage.gravitationalEnergy(sim);
    auto& totalAcousticEnergyLocal = energiesStorage.acousticEnergy(sim);
    auto& totalAcousticKineticEnergyLocal = energiesStorage.acousticKineticEnergy(sim);
    auto& totalElasticEnergyLocal = energiesStorage.elasticEnergy(sim);
    auto& totalElasticKineticEnergyLocal = energiesStorage.elasticKineticEnergy(sim);
    auto& totalPlasticMoment = energiesStorage.plasticMoment(sim);

    std::vector<Element> const& elements = meshReader->getElements();
    std::vector<Vertex> const& vertices = meshReader->getVertices();

    const auto g = seissolInstance.getGravitationSetup().acceleration;

    // Note: Default(none) is not possible, clang requires data sharing attribute for g, gcc forbids
    // it
#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for schedule(static) reduction(+ : totalGravitationalEnergyLocal,             \
                                                        totalAcousticEnergyLocal,                  \
                                                        totalAcousticKineticEnergyLocal,           \
                                                        totalElasticEnergyLocal,                   \
                                                        totalElasticKineticEnergyLocal,            \
                                                        totalMomentumX,                            \
                                                        totalMomentumY,                            \
                                                        totalMomentumZ,                            \
                                                        totalPlasticMoment)                        \
    shared(elements, vertices, lts, ltsLut, global)
#endif
    for (std::size_t elementId = 0; elementId < elements.size(); ++elementId) {
      real volume = MeshTools::volume(elements[elementId], vertices);
      CellMaterialData& material = ltsLut->lookup(lts->material, elementId);
#if defined(USE_ELASTIC) || defined(USE_VISCOELASTIC2)
      auto& cellInformation = ltsLut->lookup(lts->cellInformation, elementId);
      auto& faceDisplacements = ltsLut->lookup(lts->faceDisplacements, elementId);

      constexpr auto quadPolyDegree = CONVERGENCE_ORDER + 1;
      constexpr auto numQuadraturePointsTet = quadPolyDegree * quadPolyDegree * quadPolyDegree;

      double quadraturePointsTet[numQuadraturePointsTet][3];
      double quadratureWeightsTet[numQuadraturePointsTet];
      seissol::quadrature::TetrahedronQuadrature(
          quadraturePointsTet, quadratureWeightsTet, quadPolyDegree);

      constexpr auto numQuadraturePointsTri = quadPolyDegree * quadPolyDegree;
      double quadraturePointsTri[numQuadraturePointsTri][2];
      double quadratureWeightsTri[numQuadraturePointsTri];
      seissol::quadrature::TriangleQuadrature(
          quadraturePointsTri, quadratureWeightsTri, quadPolyDegree);

      // Needed to weight the integral.
      const auto jacobiDet = 6 * volume;

      alignas(ALIGNMENT) real numericalSolutionData[tensor::dofsQP::size()];
      auto numericalSolution = init::dofsQP::view::create(numericalSolutionData);
      // Evaluate numerical solution at quad. nodes
      kernel::evalAtQP krnl;
      krnl.evalAtQP = global->evalAtQPMatrix;
      krnl.dofsQP = numericalSolutionData;
      krnl.Q = ltsLut->lookup(lts->dofs, elementId);
      krnl.execute();

#ifdef MULTIPLE_SIMULATIONS
      auto numSub = numericalSolution.subtensor(sim, yateto::slice<>(), yateto::slice<>());
#else
      auto numSub = numericalSolution;
#endif
      for (size_t qp = 0; qp < numQuadraturePointsTet; ++qp) {
        constexpr int uIdx = 6;
        const auto curWeight = jacobiDet * quadratureWeightsTet[qp];
        const auto rho = material.local.rho;

      const auto u = numSub(qp, uIdx + 0);
      const auto v = numSub(qp, uIdx + 1);
      const auto w = numSub(qp, uIdx + 2);
      const double curKineticEnergy = 0.5 * rho * (u * u + v * v + w * w);
      const double curMomentumX = rho * u;
      const double curMomentumY = rho * v;
      const double curMomentumZ = rho * w;

      if (std::abs(material.local.mu) < 10e-14) {
        // Acoustic
        constexpr int pIdx = 0;
        const auto K = material.local.lambda;
        const auto p = numSub(qp, pIdx);
        const double curAcousticEnergy = (p * p) / (2 * K);
        totalAcousticEnergyLocal += curWeight * curAcousticEnergy;
        totalAcousticKineticEnergyLocal += curWeight * curKineticEnergy;
      } else {
        // Elastic
        totalElasticKineticEnergyLocal += curWeight * curKineticEnergy;
        auto getStressIndex = [](int i, int j) {
          const static auto lookup =
              std::array<std::array<int, 3>, 3>{{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}}};
          return lookup[i][j];
        };
        totalMomentumX += curWeight * curMomentumX;
        totalMomentumY += curWeight * curMomentumY;
        totalMomentumZ += curWeight * curMomentumZ;

        auto getStress = [&](int i, int j) { return numSub(qp, getStressIndex(i, j)); };

          const auto lambda = material.local.lambda;
          const auto mu = material.local.mu;
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
        if (cellInformation.faceTypes[face] != FaceType::freeSurfaceGravity)
          continue;

        // Displacements are stored in face-aligned coordinate system.
        // We need to rotate it to the global coordinate system.
        auto& boundaryMapping = boundaryMappings[face];
        auto Tinv = init::Tinv::view::create(boundaryMapping.TinvData);
        alignas(ALIGNMENT)
            real rotateDisplacementToFaceNormalData[init::displacementRotationMatrix::Size];

        auto rotateDisplacementToFaceNormal =
            init::displacementRotationMatrix::view::create(rotateDisplacementToFaceNormalData);
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            rotateDisplacementToFaceNormal(i, j) = Tinv(i + 6, j + 6);
          }
        }

        alignas(ALIGNMENT) std::array<real, tensor::rotatedFaceDisplacementAtQuadratureNodes::Size>
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

#ifdef MULTIPLE_SIMULATIONS
        static_assert(numQuadraturePointsTri ==
                      init::rotatedFaceDisplacementAtQuadratureNodes::Shape[1]);
        auto rotatedFaceDisplacementFused =
            init::rotatedFaceDisplacementAtQuadratureNodes::view::create(displQuadData.data());
        auto rotatedFaceDisplacement =
            rotatedFaceDisplacementFused.subtensor(sim, yateto::slice<>(), yateto::slice<>());
#else
        static_assert(numQuadraturePointsTri ==
                      init::rotatedFaceDisplacementAtQuadratureNodes::Shape[0]);
        auto rotatedFaceDisplacement =
            init::rotatedFaceDisplacementAtQuadratureNodes::view::create(displQuadData.data());
#endif
        for (unsigned i = 0; i < rotatedFaceDisplacement.shape(0); ++i) {
          // See for example (Saito, Tsunami generation and propagation, 2019) section 3.2.3 for
          // derivation.
          const auto displ = rotatedFaceDisplacement(i, 0);
          const auto curEnergy = 0.5 * rho * g * displ * displ;
          const auto curWeight = 2.0 * surface * quadratureWeightsTri[i];
          totalGravitationalEnergyLocal += curWeight * curEnergy;
        }
      }
#endif

      if (isPlasticityEnabled) {
        // plastic moment
        real* pstrainCell = ltsLut->lookup(lts->pstrain, elementId);
#ifdef USE_ANISOTROPIC
        real mu = (material.local.c44 + material.local.c55 + material.local.c66) / 3.0;
#else
        real mu = material.local.mu;
#endif
        totalPlasticMoment += mu * volume * pstrainCell[tensor::QStress::size()];
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

  const auto count = static_cast<int>(energiesStorage.energies.size());
  if (rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, energiesStorage.energies.data(), count, MPI_DOUBLE, MPI_SUM, 0, comm);
  } else {
    MPI_Reduce(energiesStorage.energies.data(),
               energiesStorage.energies.data(),
               count,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               comm);
  }
#endif
}

void EnergyOutput::reduceMinTimeSinceSlipRateBelowThreshold() {
#ifdef USE_MPI
  const auto rank = MPI::mpi.rank();
  const auto& comm = MPI::mpi.comm();

  if (rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &minTimeSinceSlipRateBelowThreshold, 1, MPI_C_REAL, MPI_MIN, 0, comm);
  } else {
    MPI_Reduce(&minTimeSinceSlipRateBelowThreshold,
               &minTimeSinceSlipRateBelowThreshold,
               1,
               MPI_C_REAL,
               MPI_MIN,
               0,
               comm);
  }
#endif
}

void EnergyOutput::printEnergies() {
  const auto rank = MPI::mpi.rank();

  if (rank == 0) {
    const auto totalAcousticEnergy =
        energiesStorage.acousticKineticEnergy() + energiesStorage.acousticEnergy();
    const auto totalElasticEnergy =
        energiesStorage.elasticKineticEnergy() + energiesStorage.elasticEnergy();
    const auto ratioElasticKinematic =
        100.0 * energiesStorage.elasticKineticEnergy() / totalElasticEnergy;
    const auto ratioElasticPotential = 100.0 * energiesStorage.elasticEnergy() / totalElasticEnergy;
    const auto ratioAcousticKinematic =
        100.0 * energiesStorage.acousticKineticEnergy() / totalAcousticEnergy;
    const auto ratioAcousticPotential =
        100.0 * energiesStorage.acousticEnergy() / totalAcousticEnergy;
    const auto totalFrictionalWork = energiesStorage.totalFrictionalWork();
    const auto staticFrictionalWork = energiesStorage.staticFrictionalWork();
    const auto radiatedEnergy = totalFrictionalWork - staticFrictionalWork;
    const auto ratioFrictionalStatic = 100.0 * staticFrictionalWork / totalFrictionalWork;
    const auto ratioFrictionalRadiated = 100.0 * radiatedEnergy / totalFrictionalWork;
    const auto ratioPlasticMoment =
        100.0 * energiesStorage.plasticMoment() /
        (energiesStorage.plasticMoment() + energiesStorage.seismicMoment());
    const auto totalMomentumX = energiesStorage.totalMomentumX();
    const auto totalMomentumY = energiesStorage.totalMomentumY();
    const auto totalMomentumZ = energiesStorage.totalMomentumZ();

      if (shouldComputeVolumeEnergies()) {
        if (totalElasticEnergy) {
          logInfo(rank) << "Elastic energy (total, % kinematic, % potential): "
                        << totalElasticEnergy << " ," << ratioElasticKinematic << " ,"
                        << ratioElasticPotential;
        }
        if (totalAcousticEnergy) {
          logInfo(rank) << "Acoustic energy (total, % kinematic, % potential): "
                        << totalAcousticEnergy << " ," << ratioAcousticKinematic << " ,"
                        << ratioAcousticPotential;
        }
        if (energiesStorage.gravitationalEnergy(sim)) {
          logInfo(rank) << "Gravitational energy:" << energiesStorage.gravitationalEnergy(sim);
        }
        if (energiesStorage.plasticMoment(sim)) {
          logInfo(rank) << "Plastic moment (value, equivalent Mw, % total moment):"
                        << energiesStorage.plasticMoment(sim) << " ,"
                        << 2.0 / 3.0 * std::log10(energiesStorage.plasticMoment(sim)) - 6.07 << " ,"
                        << ratioPlasticMoment;
        }
      } else {
        logInfo(rank) << "Volume energies skipped at this step";
      }
      if (totalAcousticEnergy) {
        logInfo(rank) << "Acoustic energy (total, % kinematic, % potential): "
                      << totalAcousticEnergy << " ," << ratioAcousticKinematic << " ,"
                      << ratioAcousticPotential;
      }
      if (energiesStorage.gravitationalEnergy()) {
        logInfo(rank) << "Gravitational energy:" << energiesStorage.gravitationalEnergy();
      }
      if (energiesStorage.plasticMoment()) {
        logInfo(rank) << "Plastic moment (value, equivalent Mw, % total moment):"
                      << energiesStorage.plasticMoment() << " ,"
                      << 2.0 / 3.0 * std::log10(energiesStorage.plasticMoment()) - 6.07 << " ,"
                      << ratioPlasticMoment;
      }
      logInfo(rank) << "Total momentum (X, Y, Z):" << totalMomentumX << " ," << totalMomentumY
                    << " ," << totalMomentumZ;
    } else {
      logInfo(rank) << "Volume energies skipped at this step";
    }

      if (totalFrictionalWork) {
        logInfo(rank) << "Frictional work (total, % static, % radiated): " << totalFrictionalWork
                      << " ," << ratioFrictionalStatic << " ," << ratioFrictionalRadiated;
        logInfo(rank) << "Seismic moment (without plasticity):"
                      << energiesStorage.seismicMoment(sim) << " Mw:"
                      << 2.0 / 3.0 * std::log10(energiesStorage.seismicMoment(sim)) - 6.07;
      }

      if (!std::isfinite(totalElasticEnergy + totalAcousticEnergy)) {
        logError() << "Detected Inf/NaN in energies. Aborting.";
      }
    }
  } }

void EnergyOutput::checkAbortCriterion(real timeSinceThreshold, const std::string& prefix_message) {
  const auto rank = MPI::mpi.rank();
  bool abort = false;
  if (rank == 0) {
    if ((timeSinceThreshold > 0) and (timeSinceThreshold < std::numeric_limits<real>::max())) {
      if (static_cast<double>(timeSinceThreshold) < terminatorMaxTimePostRupture) {
        logInfo(rank) << prefix_message.c_str() << "below threshold since" << timeSinceThreshold
                      << "s (lower than the abort criteria: " << terminatorMaxTimePostRupture
                      << "s)";
      } else {
        logInfo(rank) << prefix_message.c_str() << "below threshold since" << timeSinceThreshold
                      << "s (greater than the abort criteria: " << terminatorMaxTimePostRupture
                      << "s)";
        abort = true;
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

void EnergyOutput::writeHeader() { out << "time,variable,measurement" << std::endl; }

void EnergyOutput::writeEnergies(double time) {
  for (size_t sim = 0; sim < multipleSimulations::numberOfSimulations; sim++) {
#ifdef MULTIPLE_SIMULATIONS
    const std::string fusedSuffix = std::to_string(sim);
#else
    const std::string fusedSuffix = "";
#endif



    if (shouldComputeVolumeEnergies()) {
      out << time << ",gravitational_energy" << fusedSuffix << ","
          << energiesStorage.gravitationalEnergy(sim) << "\n"
          << time << ",acoustic_energy" << fusedSuffix << "," << energiesStorage.acousticEnergy(sim)
          << "\n"
          << time << ",acoustic_kinetic_energy" << fusedSuffix << ","
          << energiesStorage.acousticKineticEnergy(sim) << "\n"
          << time << ",elastic_energy" << fusedSuffix << "," << energiesStorage.elasticEnergy(sim)
          << "\n"
          << time << ",elastic_kinetic_energy" << fusedSuffix << ","
          << energiesStorage.elasticKineticEnergy(sim) << "\n"
          << time << ",plastic_moment" << fusedSuffix << "," << energiesStorage.plasticMoment(sim)
          << "\n"
          << time << ",momentumX," << fusedSuffix << "," << energiesStorage.totalMomentumX() << "\n"
          << time << ",momentumY," << fusedSuffix << "," << energiesStorage.totalMomentumY() << "\n"
          << time << ",momentumZ," << fusedSuffix << "," << energiesStorage.totalMomentumZ() << "\n";
    }
    out << time << ",total_frictional_work" << fusedSuffix << ","
        << energiesStorage.totalFrictionalWork(sim) << "\n"
        << time << ",static_frictional_work" << fusedSuffix << ","
        << energiesStorage.staticFrictionalWork(sim) << "\n"
        << time << ",seismic_moment" << fusedSuffix << "," << energiesStorage.seismicMoment(sim)
        << "\n"
        << time << ",potency" << fusedSuffix << "," << energiesStorage.potency(sim) << "\n"
        << time << ",plastic_moment" << fusedSuffix << "," << energiesStorage.plasticMoment(sim)
        << std::endl;
  }
}

bool EnergyOutput::shouldComputeVolumeEnergies() const {
  return outputId % computeVolumeEnergiesEveryOutput == 0;
}

} // namespace seissol::writer
