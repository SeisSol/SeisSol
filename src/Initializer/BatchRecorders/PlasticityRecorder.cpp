#include "Recorders.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void PlasticityRecorder::record(LTS &handler, Layer &layer) {
  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  setUpContext(handler, layer, loader);

  real(*pstrains)[7] = currentLayer->var(currentHandler->pstrain);
  size_t NodalStressTensorCounter = 0;
  real *scratchMem = static_cast<real *>(currentLayer->getScratchpadMemory(currentHandler->idofsScratch));
  if (currentLayer->getNumberOfCells()) {
    std::vector<real *> dofsPtrs{};
    std::vector<real *> qstressNodalPtrs{};
    std::vector<real *> pstransPtrs{};

    for (unsigned cell = 0; cell < currentLayer->getNumberOfCells(); ++cell) {
      auto Data = currentLoader->entry(cell);
      dofsPtrs.push_back(static_cast<real *>(Data.dofs));
      qstressNodalPtrs.push_back(&scratchMem[NodalStressTensorCounter]);
      NodalStressTensorCounter += tensor::QStressNodal::size();
      pstransPtrs.push_back(static_cast<real *>(pstrains[cell]));
    }
    if (!dofsPtrs.empty()) {
      ConditionalKey key(*KernelNames::Plasticity);
      checkKey(key);
      (*currentTable)[key].content[*EntityId::Dofs] = new BatchPointers(dofsPtrs);
      (*currentTable)[key].content[*EntityId::NodalStressTensor] = new BatchPointers(qstressNodalPtrs);
      (*currentTable)[key].content[*EntityId::Pstrains] = new BatchPointers(pstransPtrs);
    }
  }
}