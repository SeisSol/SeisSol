// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law::gpu {
template <typename TPMethod>
class SlipLaw : public SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod> {
  public:
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::copyStorageToLocal;

  static void copySpecificStorageDataToLocal(FrictionLawData* data,
                                             DynamicRupture::Layer& layerData) {}

  /// generic over the scalar the slip rate arrives in, so that the inversion can differentiate
  /// the state variable by the very slip rate it is solving for
  template <typename S>
  SEISSOL_DEVICE static S
      stateVariableAt(FrictionLawContext& __restrict ctx, S localSlipRate, real timeIncrement) {
    using std::exp;
    using std::expm1;
    using std::log;
    const real localSl0 = ctx.data->sl0[ctx.ltsFace][ctx.pointIndex];
    const S preexp1 = -localSlipRate * S(timeIncrement / localSl0);
    const S exp1v = exp(preexp1);
    const S exp1m = -expm1(preexp1);

    const real stateVarReference = ctx.initialVariables.stateVarReference;
    // the weighted geometric mean of L / V and Psi, with the weights 1 - e and e: an exponential
    // of logarithms, so that the large quotient and the small power never meet. As a product they
    // have to -- each factor's derivative is of the order L / V^2 against a product's of the order
    // t / V -- and in single precision the derivative then comes back as exactly zero for every
    // slip rate under a millimetre a second. The first weight goes through expm1, since
    // 1 - exp(-z) cancels wherever the relaxation is slight and the logarithm it multiplies
    // reaches eighty.
    const S logQuotient = S(std::log(localSl0)) - log(localSlipRate);
    return exp(exp1m * logQuotient + exp1v * S(std::log(stateVarReference)));
  }

  SEISSOL_DEVICE static void updateStateVariable(FrictionLawContext& __restrict ctx,
                                                 real timeIncrement) {
    ctx.stateVariableBuffer =
        stateVariableAt<real>(ctx, ctx.initialVariables.localSlipRate, timeIncrement);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_
