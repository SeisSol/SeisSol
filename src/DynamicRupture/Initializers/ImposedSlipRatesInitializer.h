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

  private:
  void rotateSlipToFaultCS(seissol::initializers::DynamicRupture* dynRup,
                           seissol::initializers::LTSTree::leaf_iterator& it,
                           real (*strikeSlip)[misc::numPaddedPoints],
                           real (*dipSlip)[misc::numPaddedPoints]);
};
} // namespace seissol::dr::initializers

#endif // SEISSOL_IMPOSEDSLIPINITIALIZER_H
