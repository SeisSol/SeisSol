#include "EnergyOutput.h"
#include <Parallel/MPI.h>
#include <Kernels/DynamicRupture.h>
#include <generated_code/kernels.h>

double seissol::writer::computeStaticWork(  GlobalData const*           global,
                                            real*                       degreesOfFreedomPlus,
                                            real*                       degreesOfFreedomMinus,
                                            DRFaceInformation const&    faceInfo,
                                            DRGodunovData const&        godunovData,
                                            real                        slip[3*seissol::model::godunovState::ld] )
{
  kernels::DynamicRupture dynRupKernel;
  
  real atQPPointsPlus[seissol::model::godunovState::reals] __attribute__((aligned(PAGESIZE_STACK)));
  real atQPPointsMinus[seissol::model::godunovState::reals] __attribute__((aligned(PAGESIZE_STACK)));

  seissol::generatedKernels::evaluateAtQuadraturePoints[4*faceInfo.plusSide](
    global->faceToNodalMatrices[faceInfo.plusSide][0],
    degreesOfFreedomPlus,
    atQPPointsPlus,
    atQPPointsPlus
  );

  seissol::generatedKernels::evaluateAtQuadraturePoints[4*faceInfo.minusSide + faceInfo.faceRelation](
    global->faceToNodalMatrices[faceInfo.minusSide][faceInfo.faceRelation],
    degreesOfFreedomMinus,
    atQPPointsMinus,
    atQPPointsMinus
  );
  
  real tractionAndSlipRate[seissol::model::godunovState::ld*seissol::model::tractionAndSlipRateMatrix::cols] __attribute__((aligned(ALIGNMENT)));
  dynRupKernel.computeTractionAndSlipRate(godunovData, atQPPointsPlus, atQPPointsMinus, tractionAndSlipRate);
  
  for (unsigned d = 0; d < 3; ++d) {
    for (unsigned i = 0; i < seissol::model::godunovState::ld; ++i) {
      tractionAndSlipRate[(3+d)*seissol::model::godunovState::ld + i] = slip[d*seissol::model::godunovState::ld + i];
    }
  }

  return -0.5 * dynRupKernel.energySpaceIntegral(godunovData, tractionAndSlipRate);
}

void seissol::writer::printEnergies(  GlobalData const*                       global,
                                      seissol::initializers::DynamicRupture*  dynRup,
                                      seissol::initializers::LTSTree*         dynRupTree)
{
  double totalWorkLocal = 0.0;
  double staticWorkLocal = 0.0;
  for (auto it = dynRupTree->beginLeaf(); it != dynRupTree->endLeaf(); ++it) {
    /// \todo timeDerivativePlus and timeDerivativeMinus are missing the last timestep.
    /// (We'd need to send the dofs over the network in order to fix this.)
    real**              timeDerivativePlus  = it->var(dynRup->timeDerivativePlus);
    real**              timeDerivativeMinus = it->var(dynRup->timeDerivativeMinus);
    DRGodunovData*      godunovData         = it->var(dynRup->godunovData);
    DRFaceInformation*  faceInformation     = it->var(dynRup->faceInformation);
    DROutput*           drOutput            = it->var(dynRup->drOutput);

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:totalWorkLocal) reduction(+:staticWorkLocal)
#endif
    for (unsigned i = 0; i < it->getNumberOfCells(); ++i) {
      if (faceInformation[i].plusSideOnThisRank) {
        totalWorkLocal += drOutput[i].frictionalEnergy;
        staticWorkLocal += computeStaticWork( global,
                                              timeDerivativePlus[i],
                                              timeDerivativeMinus[i],
                                              faceInformation[i],
                                              godunovData[i],
                                              drOutput[i].slip );
      }
    }
  }

  int rank;  
  double totalWorkGlobal = 0.0;
  double staticWorkGlobal = 0.0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI::mpi.comm(), &rank);
  MPI_Reduce(&totalWorkLocal, &totalWorkGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI::mpi.comm());
  MPI_Reduce(&staticWorkLocal, &staticWorkGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI::mpi.comm());
#else
  rank = 0;
  totalWorkGlobal = totalWorkLocal;
  staticWorkGlobal = staticWorkLocal;
#endif
  
  if (rank == 0) {
    logInfo(rank) << "Total work:" << totalWorkGlobal;
    logInfo(rank) << "Static work:" << staticWorkGlobal;
    logInfo(rank) << "Radiated energy:" << totalWorkGlobal - staticWorkGlobal;
  }
}
