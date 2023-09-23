#include "EnergyOutput.h"
#include <Common/executor.hpp>
#include <Kernels/DynamicRupture.h>
#include <Numerical_aux/Quadrature.h>
#include "DynamicRupture/Misc.h"
#include <Parallel/MPI.h>
#include "SeisSol.h"
#include "Initializer/InputParameters.hpp"

namespace seissol::writer {

double& EnergiesStorage::gravitationalEnergy() { return energies[0]; }

double& EnergiesStorage::acousticEnergy() { return energies[1]; }

double& EnergiesStorage::acousticKineticEnergy() { return energies[2]; }

double& EnergiesStorage::elasticEnergy() { return energies[3]; }

double& EnergiesStorage::elasticKineticEnergy() { return energies[4]; }

double& EnergiesStorage::totalFrictionalWork() { return energies[5]; }

double& EnergiesStorage::staticFrictionalWork() { return energies[6]; }

double& EnergiesStorage::plasticMoment() { return energies[7]; }

double& EnergiesStorage::seismicMoment() { return energies[8]; }

void EnergyOutput::init(
    seissol::initializers::GlobalDataStorage* newGlobal,
    seissol::initializers::DynRupLTSForest* dynrupForest,
    seissol::geometry::MeshReader* newMeshReader,
    seissol::initializers::ClusterLTSForest* clusterForest,
    seissol::initializers::ClusterBackmap* backmap,
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

  isFileOutputEnabled = rank == 0;
  isTerminalOutputEnabled = parameters.terminalOutput && (rank == 0);
  computeVolumeEnergiesEveryOutput = parameters.computeVolumeEnergiesEveryOutput;
  outputFileName = outputFileNamePrefix + "-energy.csv";

  global = newGlobal;
  this->dynrupForest = dynrupForest;
  meshReader = newMeshReader;
  this->clusterForest = clusterForest;
  this->backmap = backmap;

  isPlasticityEnabled = newIsPlasticityEnabled;

  Modules::registerHook(*this, SIMULATION_START);
  Modules::registerHook(*this, SYNCHRONIZATION_POINT);
  setSyncInterval(parameters.interval);
}

void EnergyOutput::syncPoint(double time) {
  assert(isEnabled);
  const auto rank = MPI::mpi.rank();
  logInfo(rank) << "Writing energy output at time" << time;
  computeEnergies();
  reduceEnergies();
  if (isTerminalOutputEnabled) {
    printEnergies();
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

template <typename Config>
static typename Config::RealT
    computeStaticWork(const typename Config::RealT* degreesOfFreedomPlus,
                      const typename Config::RealT* degreesOfFreedomMinus,
                      const DRFaceInformation& faceInfo,
                      const DRGodunovData<Config>& godunovData,
                      const typename Config::RealT
                          slip[seissol::Yateto<Config>::Tensor::slipRateInterpolated::size()],
                      const GlobalData<Config>& global) {
  using RealT = typename Config::RealT;
  RealT points[(Config::ConvergenceOrder + 1) * (Config::ConvergenceOrder + 1)][2];
  RealT spaceWeights[(Config::ConvergenceOrder + 1) * (Config::ConvergenceOrder + 1)];
  seissol::quadrature::TriangleQuadrature(points, spaceWeights, Config::ConvergenceOrder + 1);

  typename Yateto<Config>::Kernel::dynamicRupture::evaluateAndRotateQAtInterpolationPoints krnl;
  krnl.V3mTo2n = global.faceToNodalMatrices;

  alignas(PAGESIZE_STACK)
      RealT QInterpolatedPlus[Yateto<Config>::Tensor::QInterpolatedPlus::size()];
  alignas(PAGESIZE_STACK)
      RealT QInterpolatedMinus[Yateto<Config>::Tensor::QInterpolatedMinus::size()];
  alignas(Alignment)
      RealT tractionInterpolated[Yateto<Config>::Tensor::tractionInterpolated::size()];

  krnl.QInterpolated = QInterpolatedPlus;
  krnl.Q = degreesOfFreedomPlus;
  krnl.TinvT = godunovData.TinvT;
  krnl._prefetch.QInterpolated = QInterpolatedPlus;
  krnl.execute(faceInfo.plusSide, 0);

  krnl.QInterpolated = QInterpolatedMinus;
  krnl.Q = degreesOfFreedomMinus;
  krnl.TinvT = godunovData.TinvT;
  krnl._prefetch.QInterpolated = QInterpolatedMinus;
  krnl.execute(faceInfo.minusSide, faceInfo.faceRelation);

  typename Yateto<Config>::Kernel::dynamicRupture::computeTractionInterpolated trKrnl;
  trKrnl.tractionPlusMatrix = godunovData.tractionPlusMatrix;
  trKrnl.tractionMinusMatrix = godunovData.tractionMinusMatrix;
  trKrnl.QInterpolatedPlus = QInterpolatedPlus;
  trKrnl.QInterpolatedMinus = QInterpolatedMinus;
  trKrnl.tractionInterpolated = tractionInterpolated;
  trKrnl.execute();

  RealT staticFrictionalWork = 0.0;
  typename Yateto<Config>::Kernel::dynamicRupture::accumulateFrictionalEnergy feKrnl;
  feKrnl.slipRateInterpolated = slip;
  feKrnl.tractionInterpolated = tractionInterpolated;
  feKrnl.spaceWeights = spaceWeights;
  feKrnl.frictionalEnergy = &staticFrictionalWork;
  feKrnl.timeWeight = -0.5 * godunovData.doubledSurfaceArea;
  feKrnl.execute();

  return staticFrictionalWork;
}

void EnergyOutput::computeDynamicRuptureEnergies() {
  double& totalFrictionalWork = energiesStorage.totalFrictionalWork();
  double& staticFrictionalWork = energiesStorage.staticFrictionalWork();
  double& seismicMoment = energiesStorage.seismicMoment();
  dynrupForest->visitLayers([&](auto&& layerview) {
    using Config = typename std::decay_t<decltype(layerview)>::ConfigT;
    using RealT = typename Config::RealT;
    /// \todo timeDerivativePlus and timeDerivativeMinus are missing the last timestep.
    /// (We'd need to send the dofs over the network in order to fix this.)
    auto* timeDerivativePlus = layerview.var(layerview.lts.timeDerivativePlus);
    auto* timeDerivativeMinus = layerview.var(layerview.lts.timeDerivativeMinus);
    auto* godunovData = layerview.var(layerview.lts.godunovData);
    auto* faceInformation = layerview.var(layerview.lts.faceInformation);
    auto* drEnergyOutput = layerview.var(layerview.lts.drEnergyOutput);
    auto* waveSpeedsPlus = layerview.var(layerview.lts.waveSpeedsPlus);
    auto* waveSpeedsMinus = layerview.var(layerview.lts.waveSpeedsMinus);

#if defined(_OPENMP) && !defined(__NVCOMPILER)
#pragma omp parallel for reduction(                                                                \
        + : totalFrictionalWork, staticFrictionalWork, seismicMoment) default(none)                \
    shared(layerview,                                                                              \
               drEnergyOutput,                                                                     \
               faceInformation,                                                                    \
               timeDerivativeMinus,                                                                \
               timeDerivativePlus,                                                                 \
               godunovData,                                                                        \
               waveSpeedsPlus,                                                                     \
               waveSpeedsMinus)
#endif
    for (unsigned i = 0; i < layerview.layer.getNumberOfCells(); ++i) {
      if (faceInformation[i].plusSideOnThisRank) {
        for (unsigned j = 0; j < seissol::dr::misc::numberOfBoundaryGaussPoints<Config>; ++j) {
          totalFrictionalWork += drEnergyOutput[i].frictionalEnergy[j];
        }
        staticFrictionalWork += computeStaticWork(timeDerivativePlus[i],
                                                  timeDerivativeMinus[i],
                                                  faceInformation[i],
                                                  godunovData[i],
                                                  drEnergyOutput[i].slip,
                                                  global->getData<Executor::Host, Config>());

        RealT muPlus = waveSpeedsPlus[i].density * waveSpeedsPlus[i].sWaveVelocity *
                       waveSpeedsPlus[i].sWaveVelocity;
        RealT muMinus = waveSpeedsMinus[i].density * waveSpeedsMinus[i].sWaveVelocity *
                        waveSpeedsMinus[i].sWaveVelocity;
        RealT mu = 2.0 * muPlus * muMinus / (muPlus + muMinus);
        RealT seismicMomentIncrease = 0.0;
        for (unsigned k = 0; k < seissol::dr::misc::numberOfBoundaryGaussPoints<Config>; ++k) {
          seismicMomentIncrease += drEnergyOutput[i].accumulatedSlip[k];
        }
        seismicMomentIncrease *= 0.5 * godunovData[i].doubledSurfaceArea * mu /
                                 seissol::dr::misc::numberOfBoundaryGaussPoints<Config>;
        seismicMoment += seismicMomentIncrease;
      }
    }
  });
}

void EnergyOutput::computeVolumeEnergies() {
  auto& totalGravitationalEnergyLocal = energiesStorage.gravitationalEnergy();
  auto& totalAcousticEnergyLocal = energiesStorage.acousticEnergy();
  auto& totalAcousticKineticEnergyLocal = energiesStorage.acousticKineticEnergy();
  auto& totalElasticEnergyLocal = energiesStorage.elasticEnergy();
  auto& totalElasticKineticEnergyLocal = energiesStorage.elasticKineticEnergy();
  auto& totalPlasticMoment = energiesStorage.plasticMoment();

  std::vector<Element> const& elements = meshReader->getElements();
  std::vector<Vertex> const& vertices = meshReader->getVertices();

  const auto g = SeisSol::main.getGravitationSetup().acceleration;

  clusterForest->visitLayers([&](auto&& layerview) {
    using Config = typename std::decay_t<decltype(layerview)>::ConfigT;
    using MaterialT = typename Config::MaterialT;
    using RealT = typename Config::RealT;
    const auto* secondaryCellInformation = layerview.var(layerview.lts.secondaryCellInformation);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) reduction(+ : totalGravitationalEnergyLocal,             \
                                                        totalAcousticEnergyLocal,                  \
                                                        totalAcousticKineticEnergyLocal,           \
                                                        totalElasticEnergyLocal,                   \
                                                        totalElasticKineticEnergyLocal,            \
                                                        totalPlasticMoment)
#endif
    for (unsigned cell = 0; cell < layerview.layer.getNumberOfCells(); ++cell) {
      if (secondaryCellInformation[cell].duplicate == 0) {
        auto elementId = secondaryCellInformation[cell].meshId;
        RealT volume = MeshTools::volume(elements[elementId], vertices);
        if constexpr (MaterialT::Type == seissol::model::MaterialType::elastic ||
                      MaterialT::Type == seissol::model::MaterialType::viscoelastic) {
          const auto& cellInformation = layerview.var(layerview.lts.cellInformation)[cell];
          const auto& faceDisplacements = layerview.var(layerview.lts.faceDisplacements)[cell];
          const auto& material = layerview.var(layerview.lts.materialData)[cell];

          constexpr auto quadPolyDegree = Config::ConvergenceOrder + 1;
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

          alignas(Alignment) real numericalSolutionData[Yateto<Config>::Tensor::dofsQP::size()];
          auto numericalSolution =
              Yateto<Config>::Init::dofsQP::view::create(numericalSolutionData);
          // Evaluate numerical solution at quad. nodes
          typename Yateto<Config>::Kernel::evalAtQP krnl;
          krnl.evalAtQP = global->evalAtQPMatrix;
          krnl.dofsQP = numericalSolutionData;
          krnl.Q = layerview.var(layerview.lts.dofs)[cell];
          krnl.execute();

#ifdef MULTIPLE_SIMULATIONS
          auto numSub = numericalSolution.subtensor(sim, yateto::slice<>(), yateto::slice<>());
#else
          auto numSub = numericalSolution;
#endif
          for (size_t qp = 0; qp < numQuadraturePointsTet; ++qp) {
            constexpr int uIdx = 6;
            const auto curWeight = jacobiDet * quadratureWeightsTet[qp];
            const auto rho = material.rho;

            const auto u = numSub(qp, uIdx + 0);
            const auto v = numSub(qp, uIdx + 1);
            const auto w = numSub(qp, uIdx + 2);
            const double curKineticEnergy = 0.5 * rho * (u * u + v * v + w * w);

            if (std::abs(material.mu) < 10e-14) {
              // Acoustic TODO: differentiate
              constexpr int pIdx = 0;
              const auto K = material.lambda;
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
              auto getStress = [&](int i, int j) { return numSub(qp, getStressIndex(i, j)); };

              const auto lambda = material.lambda;
              const auto mu = material.mu;
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

          auto* boundaryMappings = layerview.var(layerview.lts.boundaryMapping)[cell];
          // Compute gravitational energy
          for (int face = 0; face < 4; ++face) {
            if (cellInformation.faceTypes[face] != FaceType::freeSurfaceGravity)
              continue;

            // Displacements are stored in face-aligned coordinate system.
            // We need to rotate it to the global coordinate system.
            auto& boundaryMapping = boundaryMappings[face];
            auto Tinv = Yateto<Config>::Init::Tinv::view::create(boundaryMapping.TinvData);
            alignas(Alignment) RealT rotateDisplacementToFaceNormalData
                [Yateto<Config>::Init::displacementRotationMatrix::Size];

            auto rotateDisplacementToFaceNormal =
                Yateto<Config>::Init::displacementRotationMatrix::view::create(
                    rotateDisplacementToFaceNormalData);
            for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                rotateDisplacementToFaceNormal(i, j) = Tinv(i + 6, j + 6);
              }
            }

            alignas(Alignment)
                std::array<RealT,
                           Yateto<Config>::Tensor::rotatedFaceDisplacementAtQuadratureNodes::Size>
                    displQuadData{};
            const auto* curFaceDisplacementsData = faceDisplacements[face];
            typename seissol::Yateto<
                Config>::Kernel::rotateFaceDisplacementsAndEvaluateAtQuadratureNodes evalKrnl;
            evalKrnl.rotatedFaceDisplacement = curFaceDisplacementsData;
            evalKrnl.V2nTo2JacobiQuad = Yateto<Config>::Init::V2nTo2JacobiQuad::Values;
            evalKrnl.rotatedFaceDisplacementAtQuadratureNodes = displQuadData.data();
            evalKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
            evalKrnl.execute();

            // Perform quadrature
            const auto surface = MeshTools::surface(elements[elementId], face, vertices);
            const auto rho = material.rho;

            static_assert(numQuadraturePointsTri ==
                          Yateto<Config>::Init::rotatedFaceDisplacementAtQuadratureNodes::Shape[0]);
            auto rotatedFaceDisplacement =
                Yateto<Config>::Init::rotatedFaceDisplacementAtQuadratureNodes::view::create(
                    displQuadData.data());
            for (unsigned i = 0; i < rotatedFaceDisplacement.shape(0); ++i) {
              // See for example (Saito, Tsunami generation and propagation, 2019) section 3.2.3 for
              // derivation.
              const auto displ = rotatedFaceDisplacement(i, 0);
              const auto curEnergy = 0.5 * rho * g * displ * displ;
              const auto curWeight = 2.0 * surface * quadratureWeightsTri[i];
              totalGravitationalEnergyLocal += curWeight * curEnergy;
            }
          }
        }

        if constexpr (Config::Plasticity) {
          const auto& material = layerview.var(layerview.lts.materialData)[cell];
          // plastic moment
          RealT* pstrainCell = layerview.var(layerview.lts.pstrain)[cell];
          RealT mu = material.getMu();
          totalPlasticMoment +=
              mu * volume *
              pstrainCell[6 * seissol::kernels::NumberOfAlignedBasisFunctions<RealT>(
                                  Config::ConvergenceOrder)];
        }
      }
    }
  });
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

    if (shouldComputeVolumeEnergies()) {
      if (totalElasticEnergy) {
        logInfo(rank) << "Elastic energy (total, % kinematic, % potential): " << totalElasticEnergy
                      << " ," << ratioElasticKinematic << " ," << ratioElasticPotential;
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
    } else {
      logInfo(rank) << "Volume energies skipped at this step";
    }

    if (totalFrictionalWork) {
      logInfo(rank) << "Frictional work (total, % static, % radiated): " << totalFrictionalWork
                    << " ," << ratioFrictionalStatic << " ," << ratioFrictionalRadiated;
      logInfo(rank) << "Seismic moment (without plasticity):" << energiesStorage.seismicMoment()
                    << " Mw:" << 2.0 / 3.0 * std::log10(energiesStorage.seismicMoment()) - 6.07;
    }

    if (!std::isfinite(totalElasticEnergy + totalAcousticEnergy)) {
      logError() << "Detected Inf/NaN in energies. Aborting.";
    }
  }
}

void EnergyOutput::writeHeader() { out << "time,variable,measurement" << std::endl; }

void EnergyOutput::writeEnergies(double time) {
  if (shouldComputeVolumeEnergies()) {
    out << time << ",gravitational_energy," << energiesStorage.gravitationalEnergy() << "\n"
        << time << ",acoustic_energy," << energiesStorage.acousticEnergy() << "\n"
        << time << ",acoustic_kinetic_energy," << energiesStorage.acousticKineticEnergy() << "\n"
        << time << ",elastic_energy," << energiesStorage.elasticEnergy() << "\n"
        << time << ",elastic_kinetic_energy," << energiesStorage.elasticKineticEnergy() << "\n"
        << time << ",plastic_moment," << energiesStorage.plasticMoment() << "\n";
  }
  out << time << ",total_frictional_work," << energiesStorage.totalFrictionalWork() << "\n"
      << time << ",static_frictional_work," << energiesStorage.staticFrictionalWork() << "\n"
      << time << ",seismic_moment," << energiesStorage.seismicMoment() << "\n"
      << time << ",plastic_moment," << energiesStorage.plasticMoment() << std::endl;
}

bool EnergyOutput::shouldComputeVolumeEnergies() const {
  return outputId % computeVolumeEnergiesEveryOutput == 0;
}

} // namespace seissol::writer
