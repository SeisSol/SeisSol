// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_INITIALIZER_IMPOSEDSLIPRATESINITIALIZER_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_INITIALIZER_IMPOSEDSLIPRATESINITIALIZER_H_

#include "BaseDRInitializer.h"
#include "Config.h"

namespace seissol::dr::initializer {
/**
 * Derived initializer class for the ImposedSliprates friction law
 */
class ImposedSlipRatesInitializer : public BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;

  /**
   * Fixes the parameters of the source time functions, and initializes the variables of the
   * outputs.
   */
  void deriveInitialState(DynamicRupture::Storage& drStorage) override;

  protected:
  /**
   * Reads the slip in strike and dip direction and the rupture onset, rotates the slip to the
   * coordinate systems of the faces, and sets the stresses to zero, which this friction law does
   * not use.
   */
  void readFault(DynamicRupture::Storage& drStorage) override;

  /**
   * Add additional parameters to be read from the easi file
   * This will be specialized in the derived friction law initializers
   * @param parameterToStorageMap reference to a std::unordered_map<std::string, double*>, which
   * maps the parameter name, to the address in memory, where the parameter shall be stored
   * @param layer reference to a Storage layer
   */
  void addAdditionalParameters(std::unordered_map<std::string, void*>& parameterToStorageMap,
                               DynamicRupture::Layer& layer) override = 0;

  /**
   * Ensure that all parameters are correct.
   * @param it
   */
  virtual void fixInterpolatedSTFParameters(DynamicRupture::Layer& layer);

  private:
  /**
   * Rotate slip from strike/dip cooordinate system to the fault aligned coordinate system.
   * @param it
   * @param strikeSlip: Slip in strike direction
   * @param dipSlip: Slip in dip direction
   * @param imposedSlipDirection1: Slip in fault aligned direction 1
   * @param imposedSlipDirection2: Slip in fault aligned direction 2
   */
  template <typename Cfg>
  void rotateSlipToFaultCS(
      DynamicRupture::Layer& layer,
      const std::vector<std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>>& strikeSlip,
      const std::vector<std::array<Real<Cfg>, misc::NumPaddedPoints<Cfg>>>& dipSlip,
      Real<Cfg> (*imposedSlipDirection1)[misc::NumPaddedPoints<Cfg>],
      Real<Cfg> (*imposedSlipDirection2)[misc::NumPaddedPoints<Cfg>]);
};

class ImposedSlipRatesYoffeInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;

  protected:
  void addAdditionalParameters(std::unordered_map<std::string, void*>& parameterToStorageMap,
                               DynamicRupture::Layer& layer) override;

  void fixInterpolatedSTFParameters(DynamicRupture::Layer& layer) override;
};

class ImposedSlipRatesGaussianInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;

  protected:
  void addAdditionalParameters(std::unordered_map<std::string, void*>& parameterToStorageMap,
                               DynamicRupture::Layer& layer) override;
};

class ImposedSlipRatesDeltaInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;

  protected:
  void addAdditionalParameters(std::unordered_map<std::string, void*>& parameterToStorageMap,
                               DynamicRupture::Layer& layer) override;
};

} // namespace seissol::dr::initializer

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_INITIALIZER_IMPOSEDSLIPRATESINITIALIZER_H_
