// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Vishal Sontakke (vishal.sontakke AT tum.de)
 *
 */

#ifndef SEISSOL_SRC_RESULTWRITER_POSTPROCESSOR_H_
#define SEISSOL_SRC_RESULTWRITER_POSTPROCESSOR_H_

#include "Initializer/PreProcessorMacros.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Tree/Layer.h"
#include "Initializer/Typedefs.h"
#include <array>
#include <vector>

namespace seissol {

namespace writer {

class PostProcessor {
  private:
  bool m_integrationMask[9];
  int m_numberOfVariables;
  std::vector<int> m_integerMap;
  seissol::initializer::Variable<real> m_integrals;

  public:
  PostProcessor() : m_numberOfVariables(0), m_integerMap(0L) {
    for (size_t i = 0; i < 9; i++) {
      m_integrationMask[i] = false;
    }
  }
  virtual ~PostProcessor() = default;
  void integrateQuantities(const double timestep,
                           seissol::initializer::Layer& layerData,
                           const unsigned int cell,
                           const double* const dofs);
  void setIntegrationMask(const std::array<bool, 9>& integrationMask);
  int getNumberOfVariables();
  void getIntegrationMask(bool* transferTo);
  void allocateMemory(seissol::initializer::LTSTree* ltsTree);
  const real* getIntegrals(seissol::initializer::LTSTree* ltsTree);
};

} // namespace writer

} // namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_POSTPROCESSOR_H_
