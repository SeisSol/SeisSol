// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GeneratedCode/tensor.h"
#include "Kernels/Interface.h"
#include "Recorders.h"
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <cstddef>
#include <vector>
#include <yateto.h>

using namespace device;
using namespace seissol::initializer;
using namespace seissol::initializer::recording;

void PlasticityRecorder::record(LTS::Layer& layer) {
  setUpContext(layer);

  auto* pstrains = currentLayer->var<LTS::PStrain>(AllocationPlace::Device);
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
  if (size > 0) {
    std::vector<real*> dofsPtrs(size, nullptr);
    std::vector<real*> qstressNodalPtrs(size, nullptr);
    std::vector<real*> pstransPtrs(size, nullptr);
    std::vector<real*> initialLoadPtrs(size, nullptr);
    std::vector<real*> qEtaNodalPtrs(size, nullptr);
    std::vector<real*> qStressNodalPtrs(size, nullptr);
    std::vector<real*> prevDofsPtrs(size, nullptr);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Device);
      dofsPtrs[cell] = static_cast<real*>(data.get<LTS::Dofs>());
      qstressNodalPtrs[cell] = &scratchMem[nodalStressTensorCounter];
      nodalStressTensorCounter += tensor::QStressNodal<Cfg>::size();
      pstransPtrs[cell] = static_cast<real*>(pstrains[cell]);
      initialLoadPtrs[cell] = static_cast<real*>(data.get<LTS::Plasticity>().initialLoading);
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
}
