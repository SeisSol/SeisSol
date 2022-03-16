#include "EnergyOutput.h"
#include <Kernels/DynamicRupture.h>
#include <Numerical_aux/Quadrature.h>
#include <Parallel/MPI.h>

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
  for (std::size_t meshId = 0; meshId < elements.size(); ++meshId) {
    real* pstrainCell = i_ltsLut->lookup(i_lts->pstrain, meshId);
    real volume = MeshTools::volume(elements[meshId], vertices);
    CellMaterialData& material = i_ltsLut->lookup(i_lts->material, meshId);
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

  real QInterpolatedPlus[tensor::QInterpolatedPlus::size()] alignas(PAGESIZE_STACK);
  real QInterpolatedMinus[tensor::QInterpolatedMinus::size()] alignas(PAGESIZE_STACK);
  real tractionInterpolated[tensor::tractionInterpolated::size()] alignas(ALIGNMENT);

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

  real staticWork = 0.0;
  dynamicRupture::kernel::computeFrictionalEnergy feKrnl;
  feKrnl.slipRateInterpolated = slip;
  feKrnl.tractionInterpolated = tractionInterpolated;
  feKrnl.spaceWeights = spaceWeights;
  feKrnl.frictionalEnergy = &staticWork;
  feKrnl.timeWeight = -0.5 * godunovData.doubledSurfaceArea;
  feKrnl.execute();

  return staticWork;
}

void seissol::writer::printEnergies(GlobalData const* global,
                                    seissol::initializers::DynamicRupture* dynRup,
                                    seissol::initializers::LTSTree* dynRupTree,
                                    MeshReader const& i_meshReader,
                                    seissol::initializers::LTSTree* i_ltsTree,
                                    seissol::initializers::LTS* i_lts,
                                    seissol::initializers::Lut* i_ltsLut,
                                    bool usePlasticity) {
  double totalWorkLocal = 0.0;
  double staticWorkLocal = 0.0;
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
#pragma omp parallel for reduction(+ : totalWorkLocal) reduction(+ : staticWorkLocal)
#endif
    for (unsigned i = 0; i < it->getNumberOfCells(); ++i) {
      if (faceInformation[i].plusSideOnThisRank) {
        totalWorkLocal += drOutput[i].frictionalEnergy;
        staticWorkLocal += computeStaticWork(global,
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
  double totalWorkGlobal = 0.0;
  double staticWorkGlobal = 0.0;
  double plasticMomentGlobal = 0.0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI::mpi.comm(), &rank);
  MPI_Reduce(&totalWorkLocal, &totalWorkGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI::mpi.comm());
  MPI_Reduce(&staticWorkLocal, &staticWorkGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI::mpi.comm());
  MPI_Reduce(&plasticMomentLocal, &plasticMomentGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI::mpi.comm());
#else
  rank = 0;
  totalWorkGlobal = totalWorkLocal;
  staticWorkGlobal = staticWorkLocal;
  plasticMomentGlobal = plasticMomentLocal;
#endif

  if (rank == 0) {
    logInfo(rank) << "Total work:" << totalWorkGlobal;
    logInfo(rank) << "Static work:" << staticWorkGlobal;
    logInfo(rank) << "Radiated energy:" << totalWorkGlobal - staticWorkGlobal;
    if (usePlasticity) {
      logInfo(rank) << "Total plastic moment:" << plasticMomentGlobal;
    }
  }
}
