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
  size_t nodalStressTensorCounter = 0;
  real *scratchMem = static_cast<real *>(currentLayer->getScratchpadMemory(currentHandler->idofsScratch));
  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real *> dofsPtrs(size, nullptr);
    std::vector<real *> qstressNodalPtrs(size, nullptr);
    std::vector<real *> pstransPtrs(size, nullptr);
    std::vector<real *> initialLoadPtrs(size, nullptr);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);
      dofsPtrs[cell] = static_cast<real *>(data.dofs);
      qstressNodalPtrs[cell] = &scratchMem[nodalStressTensorCounter];
      nodalStressTensorCounter += tensor::QStressNodal::size();
      pstransPtrs[cell] = static_cast<real *>(pstrains[cell]);
      initialLoadPtrs[cell] = static_cast<real *>(data.plasticity.initialLoading);
    }

    ConditionalKey key(*KernelNames::Plasticity);
    checkKey(key);
    (*currentTable)[key].content[*EntityId::Dofs] = new BatchPointers(dofsPtrs);
    (*currentTable)[key].content[*EntityId::NodalStressTensor] = new BatchPointers(qstressNodalPtrs);
    (*currentTable)[key].content[*EntityId::Pstrains] = new BatchPointers(pstransPtrs);
    (*currentTable)[key].content[*EntityId::InitialLoad] = new BatchPointers(initialLoadPtrs);
  }
}
