// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Vishal Sontakke

#include "PostProcessor.h"
#include "GeneratedCode/tensor.h"
#include <Alignment.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <array>
#include <cstddef>
#include <vector>

void seissol::writer::PostProcessor::integrateQuantities(const double timestep,
                                                         LTS::Layer& layerData,
                                                         const unsigned int cell,
                                                         const double* const dofs) {
  layerData.wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    constexpr auto NumQuantities =
        tensor::Q<Cfg>::Shape[sizeof(tensor::Q<Cfg>::Shape) / sizeof(tensor::Q<Cfg>::Shape[0]) - 1];

    // ill-defined for multisim; but irrelevant for it
    constexpr std::size_t NumAlignedBasisFunctions = tensor::Q<Cfg>::size() / NumQuantities;
    auto* integrals = layerData.var<LTS::Integrals>(cfg);
    for (int i = 0; i < m_numberOfVariables; i++) {
      integrals[cell * m_numberOfVariables + i] +=
          dofs[NumAlignedBasisFunctions * m_integerMap[i]] * timestep;
    }
  });
}

void seissol::writer::PostProcessor::setIntegrationMask(const std::vector<bool>& integrationMask) {
  m_integrationMask = integrationMask;
  for (int i = 0; i < integrationMask.size(); i++) {
    if (m_integrationMask[i]) {
      m_integerMap.push_back(i);
      m_numberOfVariables++;
    }
  }
}

int seissol::writer::PostProcessor::getNumberOfVariables() const { return m_numberOfVariables; }

void seissol::writer::PostProcessor::getIntegrationMask(bool* transferTo) {
  for (int i = 0; i < m_integrationMask.size(); i++) {
    transferTo[i] = m_integrationMask[i];
  }
}

void seissol::writer::PostProcessor::allocateMemory(LTS::Storage& ltsStorage) const {
  ltsStorage.add<LTS::Integrals>(seissol::initializer::LayerMask(Ghost),
                                 PagesizeHeap,
                                 initializer::AllocationMode::HostOnly,
                                 false,
                                 m_numberOfVariables);
}

const double* seissol::writer::PostProcessor::getIntegrals(LTS::Storage& ltsStorage) const {
  if (m_numberOfVariables == 0) {
    return nullptr;
  } else {
    return nullptr; // currently disabled: ltsStorage.var<LTS::Integrals>();
  }
}
