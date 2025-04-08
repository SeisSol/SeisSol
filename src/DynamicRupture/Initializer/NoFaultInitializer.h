// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_NOFAULTINITIALIZER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_NOFAULTINITIALIZER_H_

#include "BaseDRInitializer.h"

namespace seissol::dr::initializer {

/**
 * Derived initializer class for the NoFault friction law
 * does nothing in particular
 */
class NoFaultInitializer : public BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
  using BaseDRInitializer::initializeFault;
};
} // namespace seissol::dr::initializer

#endif // SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_NOFAULTINITIALIZER_H_
