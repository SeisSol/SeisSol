// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Vishal Sontakke

#include "PostProcessor.h"

#include "Alignment.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"

#include <array>
#include <cstddef>

void seissol::writer::PostProcessor::integrateQuantities(const double timestep,
                                                         LTS::Layer& layerData,
                                                         const unsigned int cell,
                                                         const double* const dofs) {

  real* integrals = layerData.var<LTS::Integrals>();
  constexpr auto NumQuantities =
      tensor::Q::Shape[sizeof(tensor::Q::Shape) / sizeof(tensor::Q::Shape[0]) - 1];

  // ill-defined for multisim; but irrelevant for it
  constexpr std::size_t NumAlignedBasisFunctions = tensor::Q::size() / NumQuantities;

  for (int i = 0; i < numberOfVariables_; i++) {
    integrals[cell * numberOfVariables_ + i] +=
        dofs[NumAlignedBasisFunctions * integerMap_[i]] * timestep;
  }
}

void seissol::writer::PostProcessor::setIntegrationMask(
    const std::array<bool, 9>& integrationMask) {
  for (int i = 0; i < 9; i++) {
    integrationMask_[i] = integrationMask[i];
    if (integrationMask_[i]) {
      integerMap_.push_back(i);
      numberOfVariables_++;
    }
  }
}

int seissol::writer::PostProcessor::getNumberOfVariables() const { return numberOfVariables_; }

void seissol::writer::PostProcessor::getIntegrationMask(bool* transferTo) {
  for (int i = 0; i < 9; i++) {
    transferTo[i] = integrationMask_[i];
  }
}

void seissol::writer::PostProcessor::allocateMemory(LTS::Storage& ltsStorage) const {
  ltsStorage.add<LTS::Integrals>(seissol::initializer::LayerMask(Ghost),
                                 PagesizeHeap,
                                 initializer::AllocationMode::HostOnly,
                                 false,
                                 numberOfVariables_);
}

const seissol::real* seissol::writer::PostProcessor::getIntegrals(LTS::Storage& ltsStorage) const {
  if (numberOfVariables_ == 0) {
    return nullptr;
  } else {
    return ltsStorage.var<LTS::Integrals>();
  }
}
