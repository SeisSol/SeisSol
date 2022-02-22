#include "EnergyOutput.h"
#include <Parallel/MPI.h>

void seissol::writer::printPlasticMoment(MeshReader const& i_meshReader,
                                         seissol::initializers::LTSTree* i_ltsTree,
                                         seissol::initializers::LTS* i_lts,
                                         seissol::initializers::Lut* i_ltsLut) {
  double plasticMomentLocal = 0.0;
  std::vector<Element> const& elements = i_meshReader.getElements();
  std::vector<Vertex> const& vertices = i_meshReader.getVertices();

  seissol::initializers::LayerMask ghostMask(Ghost);
  for (auto it = i_ltsTree->beginLeaf(ghostMask); it != i_ltsTree->endLeaf(); ++it) {
    real(*pstrain)[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = it->var(i_lts->pstrain);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : plasticMomentLocal)
#endif
    for (std::size_t meshId = 0; meshId < elements.size(); ++meshId) {
      real* pstrainCell = i_ltsLut->lookup(i_lts->pstrain, meshId);
      double volume = MeshTools::volume(elements[meshId], vertices);
      CellMaterialData& material = i_ltsLut->lookup(i_lts->material, meshId);
#ifdef USE_ANISOTROPIC
      double mu = (material.local.c44 + material.local.c55 + material.local.c66) / 3.0;
#else
      double mu = material.local.mu;
#endif
      plasticMomentLocal += mu * volume * pstrainCell[6 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    }
  }

  int rank;
  double plasticMomentGlobal = 0.0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI::mpi.comm(), &rank);
  MPI_Reduce(&plasticMomentLocal, &plasticMomentGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI::mpi.comm());
#else
  rank = 0;
  plasticMomentGlobal = plasticMomentLocal;
#endif

  if (rank == 0) {
    logInfo(rank) << "Total plastic moment:" << plasticMomentGlobal;
  }
}
