// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/Interface.h"
#include "Recorders.h"
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <cstddef>
#include <tensor.h>
#include <vector>
#include <yateto.h>

using namespace device;
using namespace seissol::initializer;
using namespace seissol::initializer::recording;

void PlasticityRecorder::record(LTS& handler, Layer& layer) {
  kernels::LocalData::Loader loader, loaderHost;
  loader.load(handler, layer, AllocationPlace::Device);
  loaderHost.load(handler, layer, AllocationPlace::Host);
  setUpContext(handler, layer, loader, loaderHost);

  auto* pstrains = currentLayer->var(currentHandler->pstrain, AllocationPlace::Device);
  size_t nodalStressTensorCounter = 0;
  real* scratchMem = static_cast<real*>(currentLayer->getScratchpadMemory(
      currentHandler->integratedDofsScratch, AllocationPlace::Device));
  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real*> dofsPtrs(size, nullptr);
    std::vector<real*> qstressNodalPtrs(size, nullptr);
    std::vector<real*> pstransPtrs(size, nullptr);
    std::vector<real*> initialLoadPtrs(size, nullptr);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);
      dofsPtrs[cell] = static_cast<real*>(data.dofs());
      qstressNodalPtrs[cell] = &scratchMem[nodalStressTensorCounter];
      nodalStressTensorCounter += tensor::QStressNodal::size();
      pstransPtrs[cell] = static_cast<real*>(pstrains[cell]);
      initialLoadPtrs[cell] = static_cast<real*>(data.plasticity().initialLoading);
    }

    const ConditionalKey key(*KernelNames::Plasticity);
    checkKey(key);
    (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::NodalStressTensor, qstressNodalPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::Pstrains, pstransPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::InitialLoad, initialLoadPtrs);
  }
}
