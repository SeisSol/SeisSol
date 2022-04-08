#include "EnergyOutput.h"
#include <Kernels/DynamicRupture.h>
#include <Numerical_aux/Quadrature.h>
#include <Parallel/MPI.h>
#include "SeisSol.h"

namespace seissol::writer {

double& EnergiesStorage::gravitationalEnergy() { return energies[0]; }

double& EnergiesStorage::acousticEnergy() { return energies[1]; }

double& EnergiesStorage::acousticKineticEnergy() { return energies[2]; }

double& EnergiesStorage::elasticEnergy() { return energies[3]; }

double& EnergiesStorage::elasticKineticEnergy() { return energies[4]; }

double& EnergiesStorage::totalFrictionalWork() { return energies[5]; }

double& EnergiesStorage::staticFrictionalWork() { return energies[6]; }

double& EnergiesStorage::plasticMoment() { return energies[7]; }

void EnergyOutput::init(GlobalData* newGlobal,
                        seissol::initializers::DynamicRupture* newDynRup,
                        seissol::initializers::LTSTree* newDynRuptTree,
                        MeshReader* newMeshReader,
                        seissol::initializers::LTSTree* newLtsTree,
                        seissol::initializers::LTS* newLts,
                        seissol::initializers::Lut* newLtsLut,
                        bool newIsPlasticityEnabled,
                        bool newIsTerminalOutputEnabled,
                        const std::string& outputFileNamePrefix,
                        double newSyncPointInterval) {
  if (newSyncPointInterval > 0) {
    isEnabled = true;
  } else {
    return;
  }
  const auto rank = MPI::mpi.rank();
  logInfo(rank) << "Initializing energy output.";

  isFileOutputEnabled = rank == 0;
  isTerminalOutputEnabled = newIsTerminalOutputEnabled && (rank == 0);
  outputFileName = outputFileNamePrefix + "_energy.csv";

  global = newGlobal;
  dynRup = newDynRup;
  dynRupTree = newDynRuptTree;
  meshReader = newMeshReader;
  ltsTree = newLtsTree;
  lts = newLts;
  ltsLut = newLtsLut;

  isPlasticityEnabled = newIsPlasticityEnabled;

  Modules::registerHook(*this, SIMULATION_START);
  Modules::registerHook(*this, SYNCHRONIZATION_POINT);
  setSyncInterval(newSyncPointInterval);
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
  logInfo(rank) << "Writing energy output at time" << time << "Done.";
}

void EnergyOutput::simulationStart() {
  if (isFileOutputEnabled) {
    out.open(outputFileName);
    writeHeader();
  }
  syncPoint(0.0);
}

real EnergyOutput::computePlasticMoment() {
  real plasticMoment = 0.0;
  std::vector<Element> const& elements = meshReader->getElements();
  std::vector<Vertex> const& vertices = meshReader->getVertices();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : plasticMoment) default(none) shared(elements, vertices)
#endif
  for (std::size_t elementId = 0; elementId < elements.size(); ++elementId) {
    real* pstrainCell = ltsLut->lookup(lts->pstrain, elementId);
    real volume = MeshTools::volume(elements[elementId], vertices);
    CellMaterialData& material = ltsLut->lookup(lts->material, elementId);
#ifdef USE_ANISOTROPIC
    real mu = (material.local.c44 + material.local.c55 + material.local.c66) / 3.0;
#else
    real mu = material.local.mu;
#endif
    plasticMoment += mu * volume * pstrainCell[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
  }
  return plasticMoment;
}

real EnergyOutput::computeStaticWork(const real* degreesOfFreedomPlus,
                                     const real* degreesOfFreedomMinus,
                                     const DRFaceInformation& faceInfo,
                                     const DRGodunovData& godunovData,
                                     const real slip[seissol::tensor::slipInterpolated::size()]) {
  real points[NUMBER_OF_SPACE_QUADRATURE_POINTS][2];
  real spaceWeights[NUMBER_OF_SPACE_QUADRATURE_POINTS];
  seissol::quadrature::TriangleQuadrature(points, spaceWeights, CONVERGENCE_ORDER + 1);

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints krnl;
  krnl.V3mTo2n = global->faceToNodalMatrices;

  alignas(PAGESIZE_STACK) real QInterpolatedPlus[tensor::QInterpolatedPlus::size()];
  alignas(PAGESIZE_STACK) real QInterpolatedMinus[tensor::QInterpolatedMinus::size()];
  alignas(ALIGNMENT) real tractionInterpolated[tensor::tractionInterpolated::size()];

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

  dynamicRupture::kernel::computeTractionInterpolated trKrnl;
  trKrnl.tractionPlusMatrix = godunovData.tractionPlusMatrix;
  trKrnl.tractionMinusMatrix = godunovData.tractionMinusMatrix;
  trKrnl.QInterpolatedPlus = QInterpolatedPlus;
  trKrnl.QInterpolatedMinus = QInterpolatedMinus;
  trKrnl.tractionInterpolated = tractionInterpolated;
  trKrnl.execute();

  real staticFrictionalWork = 0.0;
  dynamicRupture::kernel::accumulateFrictionalEnergy feKrnl;
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
  double& plasticMoment = energiesStorage.plasticMoment();
  for (auto it = dynRupTree->beginLeaf(); it != dynRupTree->endLeaf(); ++it) {
    /// \todo timeDerivativePlus and timeDerivativeMinus are missing the last timestep.
    /// (We'd need to send the dofs over the network in order to fix this.)
    real** timeDerivativePlus = it->var(dynRup->timeDerivativePlus);
    real** timeDerivativeMinus = it->var(dynRup->timeDerivativeMinus);
    DRGodunovData* godunovData = it->var(dynRup->godunovData);
    DRFaceInformation* faceInformation = it->var(dynRup->faceInformation);
    DROutput* drOutput = it->var(dynRup->drOutput);

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : totalFrictionalWork, staticFrictionalWork) default(none) shared(it, drOutput, faceInformation, timeDerivativeMinus, timeDerivativePlus, godunovData)
#endif
    for (unsigned i = 0; i < it->getNumberOfCells(); ++i) {
      if (faceInformation[i].plusSideOnThisRank) {
        totalFrictionalWork += drOutput[i].frictionalEnergy;
        staticFrictionalWork += computeStaticWork(timeDerivativePlus[i],
                                                  timeDerivativeMinus[i],
                                                  faceInformation[i],
                                                  godunovData[i],
                                                  drOutput[i].slip);
      }
    }
  }
  if (isPlasticityEnabled) {
    plasticMoment = computePlasticMoment();
  }
}

void EnergyOutput::computeEnergies() {
  energiesStorage.energies.fill(0.0);

  auto& totalGravitationalEnergyLocal = energiesStorage.gravitationalEnergy();
  auto& totalAcousticEnergyLocal = energiesStorage.acousticEnergy();
  auto& totalAcousticKineticEnergyLocal = energiesStorage.acousticKineticEnergy();
  auto& totalElasticEnergyLocal = energiesStorage.elasticEnergy();
  auto& totalElasticKineticEnergyLocal = energiesStorage.elasticKineticEnergy();

  std::vector<Element> const& elements = meshReader->getElements();
  std::vector<Vertex> const& vertices = meshReader->getVertices();

  int numberOfFacesReached = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) schedule(static) reduction(+ : numberOfFacesReached, totalGravitationalEnergyLocal, totalAcousticEnergyLocal, totalAcousticKineticEnergyLocal, totalElasticEnergyLocal, totalElasticKineticEnergyLocal) shared(elements, vertices, lts, ltsLut, global)
#endif
  for (std::size_t elementId = 0; elementId < elements.size(); ++elementId) {
#ifdef USE_ELASTIC
    real volume = MeshTools::volume(elements[elementId], vertices);
    CellMaterialData& material = ltsLut->lookup(lts->material, elementId);
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

      if (std::abs(material.local.mu) < 10e-14) {
        // Acoustic
        constexpr int pIdx = 0;
        const auto K = material.local.lambda;
        const auto p = numSub(qp, pIdx);

        totalAcousticEnergyLocal += 0.5 * curWeight * K * p * p;
        totalAcousticKineticEnergyLocal += curWeight * curKineticEnergy;
      } else {
        // Elastic
        totalElasticKineticEnergyLocal += curWeight * curKineticEnergy;
        auto getStressIndex = [](int i, int j) {
          auto lookup = std::array<std::array<int, 3>, 3>{{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}}};
          return lookup[i][j];
        };
        auto getStress = [&](int i, int j) { return numSub(qp, getStressIndex(i, j)); };

        const auto lambda = material.local.lambda;
        const auto mu = material.local.mu;
        const auto dilation = getStress(0, 0) + getStress(1, 1) + getStress(2, 2);
        auto computeStrain = [&](int i, int j) {
          double strain = 0.0;
          const auto factor = -1.0 * (lambda) / (2 * mu * (3 * lambda + 2 * mu));
          if (i == j) {
            strain += factor * dilation;
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
      numberOfFacesReached++;

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
      const auto g = SeisSol::main.getGravitationSetup().acceleration;

      static_assert(numQuadraturePointsTri ==
                    init::rotatedFaceDisplacementAtQuadratureNodes::Shape[0]);
      auto rotatedFaceDisplacement =
          init::rotatedFaceDisplacementAtQuadratureNodes::view::create(displQuadData.data());
      for (unsigned i = 0; i < rotatedFaceDisplacement.shape(0); ++i) {
        // See for example (Saito, Tsunami generation and propagation, 2019) section 3.2.3 for
        // derivation.
        const auto displ = rotatedFaceDisplacement(i, 0);
        const auto curEnergy = 0.5 * rho * g * displ * displ;
        const auto curWeight = 2 * surface * quadratureWeightsTri[i];
        totalGravitationalEnergyLocal += curWeight * curEnergy;
      }
    }
#endif
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
    logInfo(rank) << "Total gravitational energy:" << energiesStorage.gravitationalEnergy();
    logInfo(rank) << "Total acoustic kinetic energy:" << energiesStorage.acousticKineticEnergy();
    logInfo(rank) << "Total acoustic energy:" << energiesStorage.acousticEnergy();
    logInfo(rank) << "Total elastic kinetic energy:" << energiesStorage.elasticKineticEnergy();
    logInfo(rank) << "Total elastic strain energy:" << energiesStorage.elasticEnergy();

    const auto totalFrictionalWork = energiesStorage.totalFrictionalWork();
    const auto staticFrictionalWork = energiesStorage.staticFrictionalWork();
    const auto radiatedEnergy = totalFrictionalWork - staticFrictionalWork;
    logInfo(rank) << "Total frictional work:" << totalFrictionalWork;
    logInfo(rank) << "Static frictional work:" << staticFrictionalWork;
    logInfo(rank) << "Radiated energy:" << radiatedEnergy;
    if (isPlasticityEnabled) {
      logInfo(rank) << "Total plastic moment:" << energiesStorage.plasticMoment();
    }
  }
}

void EnergyOutput::writeHeader() {
  out << "time,"
      << "gravitational_energy,"
      << "acoustic_energy,"
      << "acoustic_kinetic_energy,"
      << "elastic_energy,"
      << "elastic_kinetic_energy,"
      << "total_frictional_work,"
      << "static_frictional_work,"
      << "plastic_moment" << std::endl;
}

void EnergyOutput::writeEnergies(double time) {
  const auto rank = MPI::mpi.rank();
  out << time << "," << energiesStorage.gravitationalEnergy() << ","
      << energiesStorage.acousticEnergy() << "," << energiesStorage.acousticKineticEnergy() << ","
      << energiesStorage.elasticEnergy() << "," << energiesStorage.elasticKineticEnergy() << ","
      << energiesStorage.totalFrictionalWork() << "," << energiesStorage.staticFrictionalWork()
      << "," << energiesStorage.plasticMoment() << std::endl;
}

} // namespace seissol::writer