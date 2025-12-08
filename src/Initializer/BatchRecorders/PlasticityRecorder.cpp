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

  layer.wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;

    auto* pstrains = currentLayer->var<LTS::PStrain>(cfg, AllocationPlace::Device);
    size_t nodalStressTensorCounter = 0;
    real* scratchMem = static_cast<real*>(
        currentLayer->var<LTS::IntegratedDofsScratch>(cfg, AllocationPlace::Device));
    real* qEtaNodalScratch =
        static_cast<real*>(currentLayer->var<LTS::QEtaNodalScratch>(cfg, AllocationPlace::Device));
    real* qStressNodalScratch = static_cast<real*>(
        currentLayer->var<LTS::QStressNodalScratch>(cfg, AllocationPlace::Device));
    real* prevDofsScratch =
        static_cast<real*>(currentLayer->var<LTS::PrevDofsScratch>(cfg, AllocationPlace::Device));
    const auto size = currentLayer->size();
    if (size > 0) {
      std::vector<void*> dofsPtrs(size, nullptr);
      std::vector<void*> qstressNodalPtrs(size, nullptr);
      std::vector<void*> pstransPtrs(size, nullptr);
      std::vector<void*> initialLoadPtrs(size, nullptr);
      std::vector<void*> qEtaNodalPtrs(size, nullptr);
      std::vector<void*> qStressNodalPtrs(size, nullptr);
      std::vector<void*> prevDofsPtrs(size, nullptr);

      for (unsigned cell = 0; cell < size; ++cell) {
        auto data = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Device);
        dofsPtrs[cell] = static_cast<real*>(data.template get<LTS::Dofs>());
        qstressNodalPtrs[cell] = &scratchMem[nodalStressTensorCounter];
        nodalStressTensorCounter += tensor::QStressNodal<Cfg>::size();
        pstransPtrs[cell] = static_cast<real*>(pstrains[cell]);
        initialLoadPtrs[cell] =
            static_cast<real*>(data.template get<LTS::Plasticity>().initialLoading);
        qEtaNodalPtrs[cell] = qEtaNodalScratch + cell * tensor::QEtaNodal<Cfg>::size();
        qStressNodalPtrs[cell] = qStressNodalScratch + cell * tensor::QStressNodal<Cfg>::size();
        prevDofsPtrs[cell] = prevDofsScratch + cell * tensor::Q<Cfg>::size();
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
  });
}
