#include "Recorders.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializer;
using namespace seissol::initializer::recording;

void PlasticityRecorder::record(LTS& handler, Layer& layer) {
  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  setUpContext(handler, layer, loader);

  auto* pstrains = currentLayer->var(currentHandler->pstrain);
  size_t nodalStressTensorCounter = 0;
  real* scratchMem =
      static_cast<real*>(currentLayer->getScratchpadMemory(currentHandler->integratedDofsScratch));
  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real*> dofsPtrs(size, nullptr);
    std::vector<real*> qstressNodalPtrs(size, nullptr);
    std::vector<real*> pstransPtrs(size, nullptr);
    std::vector<real*> initialLoadPtrs(size, nullptr);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);
      dofsPtrs[cell] = static_cast<real*>(data.dofs);
      qstressNodalPtrs[cell] = &scratchMem[nodalStressTensorCounter];
      nodalStressTensorCounter += tensor::QStressNodal::size();
      pstransPtrs[cell] = static_cast<real*>(pstrains[cell]);
      initialLoadPtrs[cell] = static_cast<real*>(data.plasticity.initialLoading);
    }

    ConditionalKey key(*KernelNames::Plasticity);
    checkKey(key);
    (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::NodalStressTensor, qstressNodalPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::Pstrains, pstransPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::InitialLoad, initialLoadPtrs);
  }
}
