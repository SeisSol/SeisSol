#include "Recorders.h"
#include "Common/configtensor.hpp"
#include <Common/cellconfigconv.hpp>
#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

template <typename Config>
void PlasticityRecorder<Config>::record(LTS<Config>& handler, Layer& layer) {
  typename kernels::LocalData<Config>::Loader loader;
  loader.load(handler, layer);
  setUpContext(handler, layer, loader);

  RealT(*pstrains)[ConfigConstants<Config>::PStrainSize] =
      currentLayer->var(currentHandler->pstrain);
  size_t nodalStressTensorCounter = 0;
  RealT* scratchMem =
      static_cast<RealT*>(currentLayer->getScratchpadMemory(currentHandler->integratedDofsScratch));
  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<RealT*> dofsPtrs(size, nullptr);
    std::vector<RealT*> qstressNodalPtrs(size, nullptr);
    std::vector<RealT*> pstransPtrs(size, nullptr);
    std::vector<RealT*> initialLoadPtrs(size, nullptr);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);
      dofsPtrs[cell] = static_cast<RealT*>(data.dofs);
      qstressNodalPtrs[cell] = &scratchMem[nodalStressTensorCounter];
      nodalStressTensorCounter += Yateto<Config>::Tensor::QStressNodal::size();
      pstransPtrs[cell] = static_cast<RealT*>(pstrains[cell]);
      initialLoadPtrs[cell] = static_cast<RealT*>(data.plasticity.initialLoading);
    }

    ConditionalKey key(*KernelNames::Plasticity);
    checkKey(key);
    (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::NodalStressTensor, qstressNodalPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::Pstrains, pstransPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::InitialLoad, initialLoadPtrs);
  }
}

const DeclareForAllConfigs<PlasticityRecorder> declPlasticityRecorder;
