// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Kernels/Interface.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Recorders.h"

#include <cstddef>
#include <vector>
#include <yateto.h>

using namespace device;
using namespace seissol::initializer;
using namespace seissol::recording;

void PlasticityRecorder::record(LTS::Layer& layer) {
  setUpContext(layer);

  size_t nodalStressTensorCounter = 0;
  real* scratchMem =
      static_cast<real*>(currentLayer->var<LTS::IntegratedDofsScratch>(AllocationPlace::Device));
  real* qEtaNodalScratch =
      static_cast<real*>(currentLayer->var<LTS::QEtaNodalScratch>(AllocationPlace::Device));
  real* qStressNodalScratch =
      static_cast<real*>(currentLayer->var<LTS::QStressNodalScratch>(AllocationPlace::Device));
  real* prevDofsScratch =
      static_cast<real*>(currentLayer->var<LTS::PrevDofsScratch>(AllocationPlace::Device));
  const auto size = currentLayer->size();

  std::size_t psize = 0;
  for (std::size_t cell = 0; cell < size; ++cell) {
    auto dataHost = currentLayer->cellRef(cell);

    if (dataHost.get<LTS::CellInformation>().plasticityEnabled) {
      ++psize;
    }
  }

  if (psize > 0) {
    std::vector<real*> dofsPtrs(psize, nullptr);
    std::vector<real*> qstressNodalPtrs(psize, nullptr);
    std::vector<real*> pstransPtrs(psize, nullptr);
    std::vector<real*> initialLoadPtrs(psize, nullptr);
    std::vector<real*> qEtaNodalPtrs(psize, nullptr);
    std::vector<real*> qStressNodalPtrs(psize, nullptr);
    std::vector<real*> prevDofsPtrs(psize, nullptr);

    std::size_t pcell = 0;
    for (std::size_t cell = 0; cell < size; ++cell) {
      const auto dataHost = currentLayer->cellRef(cell);
      auto data = currentLayer->cellRef(cell, AllocationPlace::Device);

      if (dataHost.get<LTS::CellInformation>().plasticityEnabled) {
        dofsPtrs[pcell] = static_cast<real*>(data.get<LTS::Dofs>());
        qstressNodalPtrs[pcell] = &scratchMem[nodalStressTensorCounter];
        nodalStressTensorCounter += tensor::QStressNodal::size();
        pstransPtrs[pcell] = static_cast<real*>(data.get<LTS::PStrain>());
        initialLoadPtrs[pcell] = static_cast<real*>(data.get<LTS::Plasticity>().initialLoading);
        qEtaNodalPtrs[pcell] = qEtaNodalScratch + pcell * tensor::QEtaNodal::size();
        qStressNodalPtrs[pcell] = qStressNodalScratch + pcell * tensor::QStressNodal::size();
        prevDofsPtrs[pcell] = prevDofsScratch + pcell * tensor::Q::size();
        ++pcell;
      }
    }

    const ConditionalKey key(*KernelNames::Plasticity);
    checkKey(key);
    (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::NodalStressTensor, qstressNodalPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::Pstrains, pstransPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::InitialLoad, initialLoadPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::PrevDofs, prevDofsPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::QEtaNodal, qEtaNodalPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::DuDtStrain, qStressNodalPtrs);
  }
}
