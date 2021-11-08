#ifndef SEISSOL_RATEANDSTATEINITIALIZER_H
#define SEISSOL_RATEANDSTATEINITIALIZER_H

#include "BaseDRInitializer.h"

namespace seissol::dr::initializers {
class RateAndStateInitializer;
class RateAndStateFastVelocityInitializer;
class RateAndStateThermalPressurisationInitializer;
} // namespace seissol::dr::initializers

// Aging and Slip Law share the same parameters
class seissol::dr::initializers::RateAndStateInitializer
    : public seissol::dr::initializers::BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability) override;

  protected:
  virtual void
      addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it) override;

  virtual std::pair<real, real> computeInitialStateAndFriction(real tractionXY,
                                                               real tractionXZ,
                                                               real pressure,
                                                               real rs_a,
                                                               real rs_b,
                                                               real rs_sl0,
                                                               real rs_sr0,
                                                               real rs_f0,
                                                               real initialSlipRate);
};

class seissol::dr::initializers::RateAndStateFastVelocityInitializer
    : public seissol::dr::initializers::RateAndStateInitializer {
  public:
  using RateAndStateInitializer::RateAndStateInitializer;

  protected:
  virtual void
      addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
  virtual std::pair<real, real> computeInitialStateAndFriction(real tractionXY,
                                                               real tractionXZ,
                                                               real pressure,
                                                               real rs_a,
                                                               real rs_b,
                                                               real rs_sl0,
                                                               real rs_sr0,
                                                               real rs_f0,
                                                               real initialSlipRate) override;
};

/*
 * initialize all thermal pressure parameters
 */
class seissol::dr::initializers::RateAndStateThermalPressurisationInitializer
    : public seissol::dr::initializers::RateAndStateFastVelocityInitializer {
  public:
  using RateAndStateFastVelocityInitializer::RateAndStateFastVelocityInitializer;
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree,
                               seissol::Interoperability* e_interoperability) override;

  protected:
  virtual void
      addAdditionalParameters(std::map<std::string, double*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

#endif // SEISSOL_RATEANDSTATEINITIALIZER_H
