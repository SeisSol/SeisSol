#ifndef SEISSOL_NOFAULT_H
#define SEISSOL_NOFAULT_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
  class NoFault;
}

/*
 * No friction computation
 * input stress XYStressGP, XZStressGP equals output XYTractionResultGP, XZTractionResultGP
 */
class seissol::dr::friction_law::NoFault : public BaseFrictionLaw {
public:
  virtual void evaluate(seissol::initializers::Layer &layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real fullUpdateTime,
                        real timeWeights[CONVERGENCE_ORDER]) override;
};

#endif //SEISSOL_NOFAULT_H
