// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_IMPOSEDSLIPRATESINITIALIZER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_IMPOSEDSLIPRATESINITIALIZER_H_

#include "BaseDRInitializer.h"

namespace seissol::dr::initializer {
/**
 * Derived initializer class for the ImposedSliprates friction law
 */
class ImposedSlipRatesInitializer : public BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;

  /**
   * Main function to initialize all fault dependent parameters.
   * @param drStorage pointer to the dynamic rupture storage
   * not need to store values in the Fortran parts
   */
  void initializeFault(DynamicRupture::Storage& drStorage) override;

  /**
   * Add additional parameters to be read from the easi file
   * This will be specialized in the derived friction law initializers
   * @param parameterToStorageMap reference to a std::unordered_map<std::string, double*>, which
   * maps the parameter name, to the address in memory, where the parameter shall be stored
   * @param layer reference to an Storage layer
   */
  void addAdditionalParameters(std::unordered_map<std::string, void*>& parameterToStorageMap,
                               DynamicRupture::Layer& layer) override = 0;

  /**
   * Ensure that all parameters are correct.
   * @param it
   */
  virtual void fixInterpolatedSTFParameters(DynamicRupture::Layer& layer);
};

class ImposedSlipRatesYoffeInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  void addAdditionalParameters(std::unordered_map<std::string, void*>& parameterToStorageMap,
                               DynamicRupture::Layer& layer) override;

  void fixInterpolatedSTFParameters(DynamicRupture::Layer& layer) override;
};

class ImposedSlipRatesGaussianInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  void addAdditionalParameters(std::unordered_map<std::string, void*>& parameterToStorageMap,
                               DynamicRupture::Layer& layer) override;
};

class ImposedSlipRatesDeltaInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  void addAdditionalParameters(std::unordered_map<std::string, void*>& parameterToStorageMap,
                               DynamicRupture::Layer& layer) override;
};

} // namespace seissol::dr::initializer

#endif // SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_IMPOSEDSLIPRATESINITIALIZER_H_
