// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Recorders.h"

#include <cstddef>
#include <vector>
#include <yateto.h>

using namespace seissol::initializer;
using namespace seissol::recording;

void PlasticityRecorder::record(LTS::Layer& layer) {
  setUpContext(layer);

  real* qStressNodalScratch =
      static_cast<real*>(currentLayer_->var<LTS::QStressNodalScratch>(AllocationPlace::Device));
  const auto size = currentLayer_->size();

  std::size_t psize = 0;
  for (std::size_t cell = 0; cell < size; ++cell) {
    const auto dataHost = currentLayer_->cellRef(cell);

    if (dataHost.get<LTS::CellInformation>().plasticityEnabled) {
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
      const auto dataHost = currentLayer_->cellRef(cell);
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);

      if (dataHost.get<LTS::CellInformation>().plasticityEnabled) {
        dofsPtrs[pcell] = static_cast<real*>(data.get<LTS::Dofs>());
        pstrainsPtrs[pcell] = static_cast<real*>(data.get<LTS::PStrain>());
        initialLoadPtrs[pcell] = static_cast<real*>(data.get<LTS::Plasticity>().initialLoading);
        qStressNodalPtrs[pcell] = qStressNodalScratch + pcell * tensor::QStressNodal::size();
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
