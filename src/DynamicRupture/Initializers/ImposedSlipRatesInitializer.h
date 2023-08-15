#ifndef SEISSOL_IMPOSEDSLIPINITIALIZER_H
#define SEISSOL_IMPOSEDSLIPINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
/**
 * Derived initializer class for the ImposedSliprates friction law
 */
template <typename Config>
class ImposedSlipRatesInitializer : public BaseDRInitializer<Config> {
  public:
  using BaseDRInitializer<Config>::BaseDRInitializer;

  /**
   * Main function to initialize all fault dependent parameters.
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param dynRupTree pointer to the dynamic rupture lts tree
   * not need to store values in the Fortran parts
   */
  void initializeFault(seissol::initializers::DynamicRupture<Config> const* const dynRup,
                       seissol::initializers::LTSTree* const dynRupTree) override;

  /**
   * Add additional parameters to be read from the easi file
   * This will be specialized in the derived friction law initializers
   * @param parameterToStorageMap reference to a std::unordered_map<std::string, double*>, which
   * maps the parameter name, to the address in memory, where the parameter shall be stored
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   */
  void addAdditionalParameters(
      std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
      seissol::initializers::DynamicRupture<Config> const* const dynRup,
      seissol::initializers::LTSInternalNode::leaf_iterator& it) override = 0;

  /**
   * Ensure that all parameters are correct.
   * @param dynRup
   * @param it
   */
  virtual void fixInterpolatedSTFParameters(
      seissol::initializers::DynamicRupture<Config> const* const dynRup,
      seissol::initializers::LTSInternalNode::leaf_iterator& it);

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
  void rotateSlipToFaultCS(
      seissol::initializers::DynamicRupture<Config> const* const dynRup,
      seissol::initializers::LTSTree::leaf_iterator& it,
      std::vector<std::array<typename Config::RealT, misc::numPaddedPoints<Config>>> const&
          strikeSlip,
      std::vector<std::array<typename Config::RealT, misc::numPaddedPoints<Config>>> const& dipSlip,
      typename Config::RealT (*imposedSlipDirection1)[misc::numPaddedPoints<Config>],
      typename Config::RealT (*imposedSlipDirection2)[misc::numPaddedPoints<Config>]);
};

template <typename Config>
class ImposedSlipRatesYoffeInitializer : public ImposedSlipRatesInitializer<Config> {
  using ImposedSlipRatesInitializer<Config>::ImposedSlipRatesInitializer;
  void addAdditionalParameters(
      std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
      seissol::initializers::DynamicRupture<Config> const* const dynRup,
      seissol::initializers::LTSInternalNode::leaf_iterator& it) override;

  void fixInterpolatedSTFParameters(
      seissol::initializers::DynamicRupture<Config> const* const dynRup,
      seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

template <typename Config>
class ImposedSlipRatesGaussianInitializer : public ImposedSlipRatesInitializer<Config> {
  using ImposedSlipRatesInitializer<Config>::ImposedSlipRatesInitializer;
  void addAdditionalParameters(
      std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
      seissol::initializers::DynamicRupture<Config> const* const dynRup,
      seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};
} // namespace seissol::dr::initializers

#endif // SEISSOL_IMPOSEDSLIPINITIALIZER_H
