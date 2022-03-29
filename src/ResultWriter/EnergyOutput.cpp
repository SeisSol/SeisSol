#include "EnergyOutput.h"
#include <Kernels/DynamicRupture.h>
#include <Numerical_aux/Quadrature.h>
#include <Parallel/MPI.h>
#include "SeisSol.h"

real seissol::writer::computePlasticMoment(MeshReader const& i_meshReader,
                                           seissol::initializers::LTSTree* i_ltsTree,
                                           seissol::initializers::LTS* i_lts,
                                           seissol::initializers::Lut* i_ltsLut) {
  real plasticMoment = 0.0;
  std::vector<Element> const& elements = i_meshReader.getElements();
  std::vector<Vertex> const& vertices = i_meshReader.getVertices();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : plasticMoment)
#endif
  for (std::size_t elementId = 0; elementId < elements.size(); ++elementId) {
    real* pstrainCell = i_ltsLut->lookup(i_lts->pstrain, elementId);
    real volume = MeshTools::volume(elements[elementId], vertices);
    CellMaterialData& material = i_ltsLut->lookup(i_lts->material, elementId);
#ifdef USE_ANISOTROPIC
    real mu = (material.local.c44 + material.local.c55 + material.local.c66) / 3.0;
#else
    real mu = material.local.mu;
#endif
    plasticMoment += mu * volume * pstrainCell[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
  }
  return plasticMoment;
}

real seissol::writer::computeStaticWork(GlobalData const* global,
                                        real* degreesOfFreedomPlus,
                                        real* degreesOfFreedomMinus,
                                        DRFaceInformation const& faceInfo,
                                        DRGodunovData const& godunovData,
                                        real slip[seissol::tensor::slipInterpolated::size()]) {
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

void seissol::writer::printDynamicRuptureEnergies(const GlobalData* global,
                                                  seissol::initializers::DynamicRupture* dynRup,
                                                  seissol::initializers::LTSTree* dynRupTree,
                                                  const MeshReader& i_meshReader,
                                                  seissol::initializers::LTSTree* i_ltsTree,
                                                  seissol::initializers::LTS* i_lts,
                                                  seissol::initializers::Lut* i_ltsLut, bool usePlasticity) {
  double totalFrictionalWorkLocal = 0.0;
  double staticFrictionalWorkLocal = 0.0;
  double plasticMomentLocal = 0.0;
  for (auto it = dynRupTree->beginLeaf(); it != dynRupTree->endLeaf(); ++it) {
    /// \todo timeDerivativePlus and timeDerivativeMinus are missing the last timestep.
    /// (We'd need to send the dofs over the network in order to fix this.)
    real** timeDerivativePlus = it->var(dynRup->timeDerivativePlus);
    real** timeDerivativeMinus = it->var(dynRup->timeDerivativeMinus);
    DRGodunovData* godunovData = it->var(dynRup->godunovData);
    DRFaceInformation* faceInformation = it->var(dynRup->faceInformation);
    DROutput* drOutput = it->var(dynRup->drOutput);

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : totalFrictionalWorkLocal) reduction(+ : staticFrictionalWorkLocal)
#endif
    for (unsigned i = 0; i < it->getNumberOfCells(); ++i) {
      if (faceInformation[i].plusSideOnThisRank) {
        totalFrictionalWorkLocal += drOutput[i].frictionalEnergy;
        staticFrictionalWorkLocal += computeStaticWork(global,
                                                       timeDerivativePlus[i],
                                                       timeDerivativeMinus[i],
                                                       faceInformation[i],
                                                       godunovData[i],
                                                       drOutput[i].slip);
      }
    }
  }
  if (usePlasticity) {
    plasticMomentLocal = computePlasticMoment(i_meshReader, i_ltsTree, i_lts, i_ltsLut);
  }
  int rank;
  double totalFrictionalWorkGlobal = 0.0;
  double staticFrictionalWorkGlobal = 0.0;
  double plasticMomentGlobal = 0.0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI::mpi.comm(), &rank);
  MPI_Reduce(&totalFrictionalWorkLocal,
             &totalFrictionalWorkGlobal,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             MPI::mpi.comm());
  MPI_Reduce(&staticFrictionalWorkLocal,
             &staticFrictionalWorkGlobal,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             MPI::mpi.comm());
  MPI_Reduce(&plasticMomentLocal, &plasticMomentGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI::mpi.comm());
#else
  rank = 0;
  totalFrictionalWorkGlobal = totalFrictionalWorkLocal;
  staticFrictionalWorkGlobal = staticFrictionalWorkLocal;
  plasticMomentGlobal = plasticMomentLocal;
#endif

  if (rank == 0) {
    logInfo(rank) << "Total frictional work:" << totalFrictionalWorkGlobal;
    logInfo(rank) << "Static frictional work:" << staticFrictionalWorkGlobal;
    logInfo(rank) << "Radiated energy:" << totalFrictionalWorkGlobal - staticFrictionalWorkGlobal;
    if (usePlasticity) {
      logInfo(rank) << "Total plastic moment:" << plasticMomentGlobal;
    }
  }
}

void seissol::writer::printEnergies(const GlobalData* global, seissol::initializers::DynamicRupture* dynRup,
                                    seissol::initializers::LTSTree* dynRupTree, const MeshReader& meshReader,
                                    seissol::initializers::LTSTree* ltsTree, seissol::initializers::LTS* lts,
                                    seissol::initializers::Lut* ltsLut, bool usePlasticity) {
  double totalGravitationalEnergyLocal = 0.0;
  double totalAcousticEnergyLocal = 0.0;
  double totalAcousticKineticEnergyLocal = 0.0;
  double totalElasticEnergyLocal = 0.0;
  double totalElasticKineticEnergyLocal = 0.0;

  std::vector<Element> const& elements = meshReader.getElements();
  std::vector<Vertex> const& vertices = meshReader.getVertices();

#ifdef _OPENMP
#pragma omp parallel for default(none) schedule(static) reduction(+ : totalGravitationalEnergyLocal, totalAcousticEnergyLocal, totalAcousticKineticEnergyLocal, totalElasticEnergyLocal, totalElasticKineticEnergyLocal) shared(elements, vertices, lts, ltsLut, global)
#endif
  for (std::size_t elementId = 0; elementId < elements.size(); ++elementId) {
#ifdef USE_ELASTIC
    real volume = MeshTools::volume(elements[elementId], vertices);
    CellMaterialData& material = ltsLut->lookup(lts->material, elementId);
    auto& cellInformation = ltsLut->lookup(lts->cellInformation, elementId);
    auto& faceDisplacements = ltsLut->lookup(lts->faceDisplacements, elementId);

    // Stuff below stolen from AnalysisWriter
    // TODO(Lukas) Refactor
    constexpr auto quadPolyDegree = CONVERGENCE_ORDER+1;
    constexpr auto numQuadPoints = quadPolyDegree * quadPolyDegree * quadPolyDegree;

    double quadraturePoints[numQuadPoints][3];
    double quadratureWeights[numQuadPoints];
    seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, quadPolyDegree);
    std::vector<std::array<double, 3>> quadraturePointsXyz(numQuadPoints);

    // Needed to weight the integral.
    const auto jacobiDet = 6 * volume;
    // Compute global position of quadrature points.
    double const* elementCoords[4];
    for (unsigned v = 0; v < 4; ++v) {
      elementCoords[v] = vertices[elements[elementId].vertices[ v ] ].coords;
    }
    for (unsigned int i = 0; i < numQuadPoints; ++i) {
      seissol::transformations::tetrahedronReferenceToGlobal(elementCoords[0], elementCoords[1], elementCoords[2], elementCoords[3], quadraturePoints[i], quadraturePointsXyz[i].data());
    }

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
    for (size_t qp = 0; qp < numQuadPoints; ++qp) {
      constexpr int uIdx = 6;
      const auto curWeight = jacobiDet * quadratureWeights[qp];
      const auto rho = material.local.rho;

      const auto u = numSub(qp, uIdx + 0);
      const auto v = numSub(qp, uIdx + 1);
      const auto w = numSub(qp, uIdx + 2);
      const double curKineticEnergy = 0.5 * rho * (
          u * u + v * v + w * w
      );

      if (std::abs(material.local.mu) < 10e-14 ) {
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
          auto lookup = std::array<std::array<int,3>,3>{{
              { 0, 3, 5 },
              { 3, 1, 4 },
              { 5, 4, 2 }
          }};
          return lookup[i][j];
        };
        auto getStress = [&](int i, int j) {
          return numSub(qp, getStressIndex(i,j));
        };

        const auto lambda = material.local.lambda;
        const auto mu = material.local.mu;
        const auto dilation = getStress(0,0) + getStress(1,1) + getStress(2,2);
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
            curElasticEnergy +=
                getStress(i, j) * computeStrain(i, j);
          }
        }
        totalElasticEnergyLocal += curWeight * 0.5 * curElasticEnergy;

      }
    }




      // Compute gravitational energy
    // TODO(Lukas) Fix computation of gravitational energy
    for (int face = 0; face < 4; ++face) {
      const auto surface = MeshTools::surface(elements[elementId], face, vertices);
      if (cellInformation.faceTypes[face] == FaceType::freeSurfaceGravity) {
        const auto rho = material.local.rho;
        const auto g = SeisSol::main.getGravitationSetup().acceleration;

        const auto curFaceDisplacements = faceDisplacements[face];
        // TODO(Lukas) Quadrature!
        //const auto displ = curFaceDisplacements[0];
        //totalGravitationalEnergyLocal += 0.5 * surface * rho * g * displ * displ;
      }
    }
#endif
  }


  const auto rank = MPI::mpi.rank();

  auto totalGravitationalEnergyGlobal = totalGravitationalEnergyLocal;
  //logInfo(rank) << "Total gravitational energy:" << totalGravitationalEnergyGlobal;
  logInfo(rank) << "Total acoustic kinetic energy:" << totalAcousticKineticEnergyLocal;
  logInfo(rank) << "Total acoustic energy:" << totalAcousticEnergyLocal;
  logInfo(rank) << "Total elastic kinetic energy:" << totalElasticKineticEnergyLocal;
  logInfo(rank) << "Total elastic strain energy:" << totalElasticEnergyLocal;

  return;
  printDynamicRuptureEnergies(global,
                              dynRup,
                              dynRupTree,
                              meshReader,
                              ltsTree,
                              lts,
                              ltsLut,
                              usePlasticity);
}