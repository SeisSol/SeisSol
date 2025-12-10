// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Vishal Sontakke

#ifndef SEISSOL_SRC_RESULTWRITER_POSTPROCESSOR_H_
#define SEISSOL_SRC_RESULTWRITER_POSTPROCESSOR_H_

#include "Initializer/PreProcessorMacros.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"

#include <array>
#include <vector>

namespace seissol::writer {

class PostProcessor {
  private:
  std::vector<bool> m_integrationMask;
  int m_numberOfVariables{0};
  std::vector<int> m_integerMap;

  public:
  PostProcessor() = default;
  virtual ~PostProcessor() = default;
  void integrateQuantities(double timestep,
                           LTS::Layer& layerData,
                           unsigned int cell,
                           const double* dofs);
  void setIntegrationMask(const std::vector<bool>& integrationMask);
  [[nodiscard]] int getNumberOfVariables() const;
  void getIntegrationMask(bool* transferTo);
  void allocateMemory(LTS::Storage& ltsStorage) const;
  const double* getIntegrals(LTS::Storage& ltsStorage) const;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_POSTPROCESSOR_H_
