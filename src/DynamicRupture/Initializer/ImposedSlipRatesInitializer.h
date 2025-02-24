// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
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
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param dynRupTree pointer to the dynamic rupture lts tree
   * not need to store values in the Fortran parts
   */
  void initializeFault(const seissol::initializer::DynamicRupture* dynRup,
                       seissol::initializer::LTSTree* dynRupTree) override;

  /**
   * Add additional parameters to be read from the easi file
   * This will be specialized in the derived friction law initializers
   * @param parameterToStorageMap reference to a std::unordered_map<std::string, double*>, which
   * maps the parameter name, to the address in memory, where the parameter shall be stored
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer) override = 0;

  /**
   * Ensure that all parameters are correct.
   * @param dynRup
   * @param it
   */
  virtual void fixInterpolatedSTFParameters(const seissol::initializer::DynamicRupture* dynRup,
                                            seissol::initializer::Layer& layer);

  private:
  /**
   * Rotate slip from strike/dip cooordinate system to the fault aligned coordinate system.
   * @param dynRup
   * @param it
   * @param strikeSlip: Slip in strike direction
   * @param dipSlip: Slip in dip direction
   * @param imposedSlipDirection1: Slip in fault aligned direction 1
   * @param imposedSlipDirection2: Slip in fault aligned direction 2
   */
  void rotateSlipToFaultCS(const seissol::initializer::DynamicRupture* dynRup,
                           seissol::initializer::Layer& layer,
                           const std::vector<std::array<real, misc::NumPaddedPoints>>& strikeSlip,
                           const std::vector<std::array<real, misc::NumPaddedPoints>>& dipSlip,
                           real (*imposedSlipDirection1)[misc::NumPaddedPoints],
                           real (*imposedSlipDirection2)[misc::NumPaddedPoints]);
};

class ImposedSlipRatesYoffeInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer) override;

  void fixInterpolatedSTFParameters(const seissol::initializer::DynamicRupture* dynRup,
                                    seissol::initializer::Layer& layer) override;
};

class ImposedSlipRatesGaussianInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer) override;
};

class ImposedSlipRatesDeltaInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               const seissol::initializer::DynamicRupture* dynRup,
                               seissol::initializer::Layer& layer) override;
};

} // namespace seissol::dr::initializer

#endif // SEISSOL_SRC_DYNAMICRUPTURE_INITIALIZER_IMPOSEDSLIPRATESINITIALIZER_H_
