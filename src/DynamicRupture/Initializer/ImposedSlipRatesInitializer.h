#ifndef SEISSOL_IMPOSEDSLIPINITIALIZER_H
#define SEISSOL_IMPOSEDSLIPINITIALIZER_H

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
  void initializeFault(seissol::initializer::DynamicRupture const* const dynRup,
                       seissol::initializer::LTSTree* const dynRupTree) override;

  /**
   * Add additional parameters to be read from the easi file
   * This will be specialized in the derived friction law initializers
   * @param parameterToStorageMap reference to a std::unordered_map<std::string, double*>, which
   * maps the parameter name, to the address in memory, where the parameter shall be stored
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   */
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               seissol::initializer::DynamicRupture const* const dynRup,
                               seissol::initializer::LTSInternalNode::leaf_iterator& it) override =
      0;

  /**
   * Ensure that all parameters are correct.
   * @param dynRup
   * @param it
   */
  virtual void
      fixInterpolatedSTFParameters(seissol::initializer::DynamicRupture const* const dynRup,
                                   seissol::initializer::LTSInternalNode::leaf_iterator& it);

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
  void rotateSlipToFaultCS(seissol::initializer::DynamicRupture const* const dynRup,
                           seissol::initializer::LTSTree::leaf_iterator& it,
                           std::vector<std::array<real, misc::numPaddedPoints>> const& strikeSlip,
                           std::vector<std::array<real, misc::numPaddedPoints>> const& dipSlip,
                           real (*imposedSlipDirection1)[misc::numPaddedPoints],
                           real (*imposedSlipDirection2)[misc::numPaddedPoints]);
};

class ImposedSlipRatesYoffeInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               seissol::initializer::DynamicRupture const* const dynRup,
                               seissol::initializer::LTSInternalNode::leaf_iterator& it) override;

  void fixInterpolatedSTFParameters(
      seissol::initializer::DynamicRupture const* const dynRup,
      seissol::initializer::LTSInternalNode::leaf_iterator& it) override;
};

class ImposedSlipRatesGaussianInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  void addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                               seissol::initializer::DynamicRupture const* const dynRup,
                               seissol::initializer::LTSInternalNode::leaf_iterator& it) override;
};
} // namespace seissol::dr::initializer

#endif // SEISSOL_IMPOSEDSLIPINITIALIZER_H
