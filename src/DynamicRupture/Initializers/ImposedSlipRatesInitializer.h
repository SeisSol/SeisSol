#ifndef SEISSOL_IMPOSEDSLIPINITIALIZER_H
#define SEISSOL_IMPOSEDSLIPINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
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
   * @param e_interoperability pointer to the interoperability instance, can be removed once we do
   * not need to store values in the Fortran parts
   */
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* eInteroperability);

  /**
   * Add additional parameters to be read from the easi file
   * This will be specialized in the derived friction law initializers
   * @param parameterToStorageMap reference to a std::unordered_map<std::string, double*>, which
   * maps the parameter name, to the address in memory, where the parameter shall be stored
   * @param dynRup pointer to the respective dynamic rupture datastructure
   * @param it reference to an LTSTree leaf_iterator
   */
  virtual void
      addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it) = 0;

  /**
   * Ensure that all parameters are correct.
   * @param dynRup
   * @param it
   */
  virtual void
      fixInterpolatedSTFParameters(seissol::initializers::DynamicRupture* dynRup,
                                   seissol::initializers::LTSInternalNode::leaf_iterator& it);

  private:
  /**
   * Rotate slip from strike/dip cooordinate system to the fault aligned coordinate system.
   * @param dynRup
   * @param it
   * @param strikeSlip: Slip in strike direction
   * @param dipSlip: Slip in dip direction
   * @param imposedSlipDirections: Slip in fault aligned directions
   */
  void rotateSlipToFaultCS(seissol::initializers::DynamicRupture* dynRup,
                           seissol::initializers::LTSTree::leaf_iterator& it,
                           std::vector<std::array<real, misc::numPaddedPoints>> const& strikeSlip,
                           std::vector<std::array<real, misc::numPaddedPoints>> const& dipSlip,
                           real (*imposedSlipDirections)[misc::numPaddedPoints][2]);
};

class ImposedSlipRatesYoffeInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  virtual void
      addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it);

  virtual void fixInterpolatedSTFParameters(
      seissol::initializers::DynamicRupture* dynRup,
      seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

class ImposedSlipRatesGaussianInitializer : public ImposedSlipRatesInitializer {
  using ImposedSlipRatesInitializer::ImposedSlipRatesInitializer;
  virtual void
      addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it);
};
} // namespace seissol::dr::initializers

#endif // SEISSOL_IMPOSEDSLIPINITIALIZER_H
