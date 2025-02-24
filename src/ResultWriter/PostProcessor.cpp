// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Vishal Sontakke

#include "PostProcessor.h"
#include <Common/Constants.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Layer.h>
#include <array>

void seissol::writer::PostProcessor::integrateQuantities(const double timestep,
                                                         seissol::initializer::Layer& layerData,
                                                         const unsigned int cell,
                                                         const double* const dofs) {

  real* integrals = layerData.var(m_integrals);
  for (int i = 0; i < m_numberOfVariables; i++) {
    integrals[cell * m_numberOfVariables + i] +=
        dofs[NumAlignedBasisFunctions * m_integerMap[i]] * timestep;
  }
}

void seissol::writer::PostProcessor::setIntegrationMask(
    const std::array<bool, 9>& integrationMask) {
  for (int i = 0; i < 9; i++) {
    m_integrationMask[i] = integrationMask[i];
    if (m_integrationMask[i]) {
      m_integerMap.push_back(i);
      m_numberOfVariables++;
    }
  }
  m_integrals.count = m_numberOfVariables;
}

int seissol::writer::PostProcessor::getNumberOfVariables() const { return m_numberOfVariables; }

void seissol::writer::PostProcessor::getIntegrationMask(bool* transferTo) {
  for (int i = 0; i < 9; i++) {
    transferTo[i] = m_integrationMask[i];
  }
}

void seissol::writer::PostProcessor::allocateMemory(seissol::initializer::LTSTree* ltsTree) {
  ltsTree->addVar(m_integrals,
                  seissol::initializer::LayerMask(Ghost),
                  PagesizeHeap,
                  initializer::AllocationMode::HostOnly);
}

const real* seissol::writer::PostProcessor::getIntegrals(seissol::initializer::LTSTree* ltsTree) {
  if (m_numberOfVariables == 0) {
    return nullptr;
  } else {
    return ltsTree->var(m_integrals);
  }
}
