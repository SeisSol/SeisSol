// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_LINEARSLIPWEAKENINGINITIALIZER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_LINEARSLIPWEAKENINGINITIALIZER_H_

#include "BaseDRInitializer.h"

namespace seissol::dr::initializer {

/**
 * Derived initializer class for the LinearSlipWeakening friction law
 */
class LinearSlipWeakeningInitializer : public BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
  /**
   * Computes initial friction and slip rates
   */
  void initializeFault(const seissol::initializer::DynamicRupture* dynRup,
                       seissol::initializer::LTSTree* dynRupTree) override;

  protected:
  /**
   * Adds the additional parameters mu_s, mu_d, d_c, cohesion and if available forced_rupture_time.
   */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer) override;
};

/**
 * Derived initializer class for the LinearSlipWeakening friction law with bimaterial regularization
 */
class LinearSlipWeakeningBimaterialInitializer : public LinearSlipWeakeningInitializer {
  public:
  using LinearSlipWeakeningInitializer::LinearSlipWeakeningInitializer;
  /**
   * Computes initial value for the regularized strength
   */
  void initializeFault(const seissol::initializer::DynamicRupture* dynRup,
                       seissol::initializer::LTSTree* dynRupTree) override;
};

} // namespace seissol::dr::initializer

#endif // SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_LINEARSLIPWEAKENINGINITIALIZER_H_
