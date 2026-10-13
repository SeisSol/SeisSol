// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Recorders.h"

#include <cstddef>
#include <vector>
#include <yateto.h>

using namespace seissol::initializer;
using namespace seissol::recording;

template <typename Cfg>
void PlasticityRecorder<Cfg>::record(LTS::Layer& layer) {
  setUpContext(layer);

  real* qStressNodalScratch = static_cast<real*>(
      currentLayer_->var<LTS::QStressNodalScratch>(Cfg(), AllocationPlace::Device));
  const auto size = currentLayer_->size();

  std::size_t psize = 0;
  for (std::size_t cell = 0; cell < size; ++cell) {
    const auto dataHost = currentLayer_->cellRef<Cfg>(cell);

    if (dataHost.template get<LTS::CellInformation>().plasticityEnabled) {
      ++psize;
    }
  }

  if (psize > 0) {
    std::vector<real*> dofsPtrs(psize, nullptr);
    std::vector<real*> pstrainsPtrs(psize, nullptr);
    std::vector<real*> initialLoadPtrs(psize, nullptr);
    std::vector<real*> qStressNodalPtrs(psize, nullptr);

    std::size_t pcell = 0;
    for (std::size_t cell = 0; cell < size; ++cell) {
      const auto dataHost = currentLayer_->cellRef<Cfg>(cell);
      auto data = currentLayer_->cellRef<Cfg>(cell, AllocationPlace::Device);

      if (dataHost.template get<LTS::CellInformation>().plasticityEnabled) {
        dofsPtrs[pcell] = static_cast<real*>(data.template get<LTS::Dofs>());
        pstrainsPtrs[pcell] = static_cast<real*>(data.template get<LTS::PStrain>());
        initialLoadPtrs[pcell] =
            static_cast<real*>(data.template get<LTS::Plasticity>().initialLoading);
        qStressNodalPtrs[pcell] = qStressNodalScratch + pcell * tensor::QStressNodal<Cfg>::size();
        ++pcell;
      }
    }

    const ConditionalKey key(*KernelNames::Plasticity);
    checkKey(key);
    (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::NodalStressTensor, qStressNodalPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::Pstrains, pstrainsPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::InitialLoad, initialLoadPtrs);
  }
}

#define SEISSOL_INSTANTIATE(Cfg) template class seissol::recording::PlasticityRecorder<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE
