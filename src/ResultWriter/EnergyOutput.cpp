#include "EnergyOutput.h"
#include <Kernels/DynamicRupture.h>
#include <Numerical_aux/Quadrature.h>
#include <Parallel/MPI.h>
#include "SeisSol.h"
#include <sstream>

namespace seissol::writer {

double& EnergiesStorage::gravitationalEnergy(size_t sim) {
  return energies[0 + sim * 9];
}
double& EnergiesStorage::acousticEnergy(size_t sim) {
  return energies[1 + sim * 9];
}
double& EnergiesStorage::acousticKineticEnergy(size_t sim) {
  return energies[2 + sim * 9];
}
double& EnergiesStorage::elasticEnergy(size_t sim) {
  return energies[3 + sim * 9];
}
double& EnergiesStorage::elasticKineticEnergy(size_t sim) {
  return energies[4 + sim * 9];
}
double& EnergiesStorage::totalFrictionalWork(size_t sim) {
  return energies[5 + sim * 9];
}
double& EnergiesStorage::staticFrictionalWork(size_t sim) {
  return energies[6 + sim * 9];
}
double& EnergiesStorage::plasticMoment(size_t sim) {
  return energies[7 + sim * 9];
}
double& EnergiesStorage::seismicMoment(size_t sim) {
  return energies[8 + sim * 9];
}

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
  outputFileName = outputFileNamePrefix + "-energy.csv";

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

std::array<real, multipleSimulations::numberOfSimulations> EnergyOutput::computeStaticWork(const real* degreesOfFreedomPlus,
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

  alignas(ALIGNMENT) real staticFrictionalWork[tensor::frictionalEnergy::size()]{};

  dynamicRupture::kernel::accumulateFrictionalEnergy feKrnl;
  feKrnl.slipRateInterpolated = slip;
  feKrnl.tractionInterpolated = tractionInterpolated;
  feKrnl.spaceWeights = spaceWeights;
  feKrnl.frictionalEnergy = staticFrictionalWork;
  feKrnl.timeWeight = -0.5 * godunovData.doubledSurfaceArea;
  feKrnl.execute();

  std::array<real, multipleSimulations::numberOfSimulations> frictionalWorkReturn;
  std::copy(staticFrictionalWork, staticFrictionalWork + multipleSimulations::numberOfSimulations, std::begin(frictionalWorkReturn));
  return frictionalWorkReturn;
}

void EnergyOutput::computeDynamicRuptureEnergies() {
  for (size_t s = 0; s < multipleSimulations::numberOfSimulations; s++) {
    double& totalFrictionalWork= energiesStorage.totalFrictionalWork(s);
    double& staticFrictionalWork= energiesStorage.staticFrictionalWork(s);
    double& seismicMoment= energiesStorage.seismicMoment(s);

    for (auto it = dynRupTree->beginLeaf(); it != dynRupTree->endLeaf(); ++it) {
      /// \todo timeDerivativePlus and timeDerivativeMinus are missing the last timestep.
      /// (We'd need to send the dofs over the network in order to fix this.)
      real** timeDerivativePlus = it->var(dynRup->timeDerivativePlus);
      real** timeDerivativeMinus = it->var(dynRup->timeDerivativeMinus);
      DRGodunovData* godunovData = it->var(dynRup->godunovData);
      DRFaceInformation* faceInformation = it->var(dynRup->faceInformation);
      DROutput* drOutput = it->var(dynRup->drOutput);
      seissol::model::IsotropicWaveSpeeds* waveSpeedsPlus = it->var(dynRup->waveSpeedsPlus);
      seissol::model::IsotropicWaveSpeeds* waveSpeedsMinus = it->var(dynRup->waveSpeedsMinus);

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : totalFrictionalWork, staticFrictionalWork, seismicMoment) default(none) shared(s, it, drOutput, faceInformation, timeDerivativeMinus, timeDerivativePlus, godunovData, waveSpeedsPlus, waveSpeedsMinus)
#endif
      for (unsigned i = 0; i < it->getNumberOfCells(); ++i) {
        if (faceInformation[i].plusSideOnThisRank) {
          totalFrictionalWork += drOutput[i].frictionalEnergy;
          staticFrictionalWork += computeStaticWork(timeDerivativePlus[i],
                                                    timeDerivativeMinus[i],
                                                    faceInformation[i],
                                                    godunovData[i],
                                                    drOutput[i].slip)[s];

          real muPlus = waveSpeedsPlus[i].density * waveSpeedsPlus[i].sWaveVelocity *
                        waveSpeedsPlus[i].sWaveVelocity;
          real muMinus = waveSpeedsMinus[i].density * waveSpeedsMinus[i].sWaveVelocity *
                         waveSpeedsMinus[i].sWaveVelocity;
          real mu = 2.0 * muPlus * muMinus / (muPlus + muMinus);
          real seismicMomentIncrease = 0.0;
          for (unsigned k = 0; k < seissol::tensor::squaredNormSlipRateInterpolated::Shape[0]; ++k) {
            seismicMomentIncrease += drOutput[i].accumulatedSlip[k];
          }
          seismicMomentIncrease *= 0.5 * godunovData[i].doubledSurfaceArea * mu /
                                   seissol::tensor::squaredNormSlipRateInterpolated::Shape[0];
          seismicMoment += seismicMomentIncrease;
        }
      }
    }
  }
}

void EnergyOutput::computeEnergies() {
  energiesStorage.energies.fill(0.0);

  for (size_t s = 0; s < multipleSimulations::numberOfSimulations; s++) {
    auto& totalGravitationalEnergyLocal = energiesStorage.gravitationalEnergy(s);
    auto& totalAcousticEnergyLocal = energiesStorage.acousticEnergy(s);
    auto& totalAcousticKineticEnergyLocal = energiesStorage.acousticKineticEnergy(s);
    auto& totalElasticEnergyLocal = energiesStorage.elasticEnergy(s);
    auto& totalElasticKineticEnergyLocal = energiesStorage.elasticKineticEnergy(s);
    auto& totalPlasticMoment = energiesStorage.plasticMoment(s);
  
    std::vector<Element> const& elements = meshReader->getElements();
    std::vector<Vertex> const& vertices = meshReader->getVertices();
  
    const auto g = SeisSol::main.getGravitationSetup().acceleration;
  
    // Note: Default(none) is not possible, clang requires data sharing attribute for g, gcc forbids
    // it
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : totalGravitationalEnergyLocal, totalAcousticEnergyLocal, totalAcousticKineticEnergyLocal, totalElasticEnergyLocal, totalElasticKineticEnergyLocal, totalPlasticMoment) shared(elements, vertices, lts, ltsLut, global)
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
      auto numSub = numericalSolution.subtensor(s, yateto::slice<>(), yateto::slice<>());
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

          const double curAcousticEnergy = (p * p) / (2 * K);
          totalAcousticEnergyLocal += curWeight * curAcousticEnergy;
          totalAcousticKineticEnergyLocal += curWeight * curKineticEnergy;
        } else {
          // Elastic
          totalElasticKineticEnergyLocal += curWeight * curKineticEnergy;
          auto getStressIndex = [](int i, int j) {
            const static auto lookup = std::array<std::array<int, 3>, 3>{{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}}};
            return lookup[i][j];
          };
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

        static_assert(numQuadraturePointsTri ==
                      init::rotatedFaceDisplacementAtQuadratureNodes::Shape[multipleSimulations::basisFunctionDimension]);

        auto rotatedFaceDisplacement =
            init::rotatedFaceDisplacementAtQuadratureNodes::view::create(displQuadData.data());
#ifdef MULTIPLE_SIMULATIONS
        auto rotatedSub = rotatedFaceDisplacement.subtensor(s, yateto::slice<>(), yateto::slice<>());
#else
        auto rotatedSub = rotatedFaceDisplacement;
#endif
        for (unsigned i = 0; i < rotatedFaceDisplacement.shape(0); ++i) {
          // See for example (Saito, Tsunami generation and propagation, 2019) section 3.2.3 for
          // derivation.
          const auto displ = rotatedSub(i, 0);
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
        totalPlasticMoment += mu * volume * pstrainCell[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
      }
    }
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
    for (size_t s = 0; s < multipleSimulations::numberOfSimulations; s++) {
      const auto totalAcousticEnergy =
          energiesStorage.acousticKineticEnergy(s) + energiesStorage.acousticEnergy(s);
      const auto totalElasticEnergy =
          energiesStorage.elasticKineticEnergy(s) + energiesStorage.elasticEnergy(s);
      const auto ratioElasticKinematic =
          100.0 * energiesStorage.elasticKineticEnergy(s) / totalElasticEnergy;
      const auto ratioElasticPotential = 100.0 * energiesStorage.elasticEnergy(s) / totalElasticEnergy;
      const auto ratioAcousticKinematic =
          100.0 * energiesStorage.acousticKineticEnergy(s) / totalAcousticEnergy;
      const auto ratioAcousticPotential =
          100.0 * energiesStorage.acousticEnergy(s) / totalAcousticEnergy;
      const auto totalFrictionalWork = energiesStorage.totalFrictionalWork();
      const auto staticFrictionalWork = energiesStorage.staticFrictionalWork();
      const auto radiatedEnergy = totalFrictionalWork - staticFrictionalWork;
      const auto ratioFrictionalStatic = 100.0 * staticFrictionalWork / totalFrictionalWork;
      const auto ratioFrictionalRadiated = 100.0 * radiatedEnergy / totalFrictionalWork;
      const auto ratioPlasticMoment =
          100.0 * energiesStorage.plasticMoment(s) /
          (energiesStorage.plasticMoment(s) + energiesStorage.seismicMoment(s));
#ifdef MULTIPLE_SIMULATIONS
      std::stringstream multipleSimulationStream;
      multipleSimulationStream << " for simulation " << s;
      const std::string multipleSimulationString = multipleSimulationStream.str();
      const char* multipleSimulationSnippet = multipleSimulationString.c_str();
#else
      //note: utils/logger puts std::string between quotation marks.
      const char* multipleSimulationSnippet = "";
#endif
      if (totalElasticEnergy) {
        logInfo(rank) << "Elastic energy" << multipleSimulationSnippet
                      << "(total, % kinematic, % potential): " << totalElasticEnergy
                      << " ," << ratioElasticKinematic << " ," << ratioElasticPotential;
      }
      if (totalAcousticEnergy) {
        logInfo(rank) << "Acoustic energy" << multipleSimulationSnippet
                      << "(total, % kinematic, % potential): " << totalAcousticEnergy
                      << " ," << ratioAcousticKinematic << " ," << ratioAcousticPotential;
      }
      if (energiesStorage.gravitationalEnergy(s)) {
        logInfo(rank) << "Gravitational energy" << multipleSimulationSnippet
                      << ":" << energiesStorage.gravitationalEnergy(s);
      }
      if (totalFrictionalWork) {
        logInfo(rank) << "Frictional work" << multipleSimulationSnippet
                      << "(total, % static, % radiated): " << totalFrictionalWork
                      << " ," << ratioFrictionalStatic << " ," << ratioFrictionalRadiated;
        logInfo(rank) << "Seismic moment" << multipleSimulationSnippet
                      << "(without plasticity):" << energiesStorage.seismicMoment()
                      << " Mw:" << 2.0 / 3.0 * std::log10(energiesStorage.seismicMoment()) - 6.07;
      }
      if (energiesStorage.plasticMoment(s)) {
        logInfo(rank) << "Plastic moment" << multipleSimulationSnippet
                      << "(value, equivalent Mw, % total moment):"
                      << energiesStorage.plasticMoment(s) << " ,"
                      << 2.0 / 3.0 * std::log10(energiesStorage.plasticMoment()) - 6.07 << " ,"
                      << ratioPlasticMoment;
        ;
      }
    }
  }
}

void EnergyOutput::writeHeader() {
  out << "time,"
#ifdef MULTIPLE_SIMULATIONS
      << "simulation_index,"
#endif
      << "gravitational_energy,"
      << "acoustic_energy,"
      << "acoustic_kinetic_energy,"
      << "elastic_energy,"
      << "elastic_kinetic_energy,"
      << "total_frictional_work,"
      << "static_frictional_work,"
      << "seismic_moment,"
      << "plastic_moment" << std::endl;
}

void EnergyOutput::writeEnergies(double time) {
  for (size_t s = 0; s < multipleSimulations::numberOfSimulations; s++) {
    out << time
#ifdef MULTIPLE_SIMULATIONS
        << "," << s
#endif
        << "," << energiesStorage.gravitationalEnergy(s) << ","
        << energiesStorage.acousticEnergy(s) << "," << energiesStorage.acousticKineticEnergy(s) << ","
        << energiesStorage.elasticEnergy(s) << "," << energiesStorage.elasticKineticEnergy(s) << ","
        << energiesStorage.totalFrictionalWork(s) << "," << energiesStorage.staticFrictionalWork(s)
        << "," << energiesStorage.seismicMoment(s) << "," << energiesStorage.plasticMoment(s)
        << std::endl;
  }
}

} // namespace seissol::writer
