#include "RateAndStateInitializer.h"

#include "DynamicRupture/Misc.h"

namespace seissol::dr::initializers {
template <typename Config>
void RateAndStateInitializer<Config>::initializeFault(
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  BaseDRInitializer<Config>::initializeFault(dynRup, dynRupTree);
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSRateAndState<Config> const* const>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    bool(*dynStressTimePending)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->dynStressTimePending);
    typename Config::RealT(*slipRate1)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->slipRate1);
    typename Config::RealT(*slipRate2)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->slipRate2);
    typename Config::RealT(*mu)[misc::numPaddedPoints<Config>] = it->var(concreteLts->mu);

    typename Config::RealT(*stateVariable)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->stateVariable);
    typename Config::RealT(*rsSl0)[misc::numPaddedPoints<Config>] = it->var(concreteLts->rsSl0);
    typename Config::RealT(*rsA)[misc::numPaddedPoints<Config>] = it->var(concreteLts->rsA);
    typename Config::RealT(*initialStressInFaultCS)[misc::numPaddedPoints<Config>][6] =
        it->var(concreteLts->initialStressInFaultCS);

    const typename Config::RealT initialSlipRate =
        misc::magnitude(drParameters->rsInitialSlipRate1, drParameters->rsInitialSlipRate2);

    using namespace dr::misc::quantity_indices;
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = drParameters->isDsOutputOn;
        slipRate1[ltsFace][pointIndex] = drParameters->rsInitialSlipRate1;
        slipRate2[ltsFace][pointIndex] = drParameters->rsInitialSlipRate2;
        // compute initial friction and state
        auto stateAndFriction =
            computeInitialStateAndFriction(initialStressInFaultCS[ltsFace][pointIndex][XY],
                                           initialStressInFaultCS[ltsFace][pointIndex][XZ],
                                           initialStressInFaultCS[ltsFace][pointIndex][XX],
                                           rsA[ltsFace][pointIndex],
                                           drParameters->rsB,
                                           rsSl0[ltsFace][pointIndex],
                                           drParameters->rsSr0,
                                           drParameters->rsF0,
                                           initialSlipRate);
        stateVariable[ltsFace][pointIndex] = stateAndFriction.stateVariable;
        mu[ltsFace][pointIndex] = stateAndFriction.frictionCoefficient;
      }
    }
  }
}

template <typename Config>
typename RateAndStateInitializer<Config>::StateAndFriction
    RateAndStateInitializer<Config>::computeInitialStateAndFriction(
        typename Config::RealT traction1,
        typename Config::RealT traction2,
        typename Config::RealT pressure,
        typename Config::RealT rsA,
        typename Config::RealT rsB,
        typename Config::RealT rsSl0,
        typename Config::RealT rsSr0,
        typename Config::RealT rsF0,
        typename Config::RealT initialSlipRate) {
  StateAndFriction result;
  const double absoluteTraction = misc::magnitude(traction1, traction2);
  const double tmp = std::abs(absoluteTraction / (rsA * pressure));
  result.stateVariable = rsSl0 / rsSr0 *
                         std::exp((rsA * std::log(std::exp(tmp) - std::exp(-tmp)) - rsF0 -
                                   rsA * std::log(initialSlipRate / rsSr0)) /
                                  rsB);
  if (result.stateVariable < 0) {
    logWarning()
        << "Found a negative state variable while initializing the fault. Are you sure your "
           "setup is correct?";
  }
  const double tmp2 = initialSlipRate * 0.5 / rsSr0 *
                      std::exp((rsF0 + rsB * std::log(rsSr0 * result.stateVariable / rsSl0)) / rsA);
  result.frictionCoefficient = rsA * std::asinh(tmp2);
  return result;
}

template <typename Config>
void RateAndStateInitializer<Config>::addAdditionalParameters(
    std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSRateAndState<Config> const* const>(dynRup);
  typename Config::RealT(*rsSl0)[misc::numPaddedPoints<Config>] = it->var(concreteLts->rsSl0);
  typename Config::RealT(*rsA)[misc::numPaddedPoints<Config>] = it->var(concreteLts->rsA);
  parameterToStorageMap.insert({"rs_sl0", (typename Config::RealT*)rsSl0});
  parameterToStorageMap.insert({"rs_a", (typename Config::RealT*)rsA});
}

template <typename Config>
typename RateAndStateInitializer<Config>::StateAndFriction
    RateAndStateFastVelocityInitializer<Config>::computeInitialStateAndFriction(
        typename Config::RealT traction1,
        typename Config::RealT traction2,
        typename Config::RealT pressure,
        typename Config::RealT rsA,
        typename Config::RealT rsB,
        typename Config::RealT rsSl0,
        typename Config::RealT rsSr0,
        typename Config::RealT rsF0,
        typename Config::RealT initialSlipRate) {
  typename RateAndStateInitializer<Config>::StateAndFriction result;
  const typename Config::RealT absoluteTraction = misc::magnitude(traction1, traction2);
  const typename Config::RealT tmp = std::abs(absoluteTraction / (rsA * pressure));
  result.stateVariable =
      rsA * std::log(2.0 * rsSr0 / initialSlipRate * (std::exp(tmp) - std::exp(-tmp)) / 2.0);
  if (result.stateVariable < 0) {
    logWarning()
        << "Found a negative state variable while initializing the fault. Are you sure your "
           "setup is correct?";
  }
  const typename Config::RealT tmp2 =
      initialSlipRate * 0.5 / rsSr0 * std::exp(result.stateVariable / rsA);
  result.frictionCoefficient = rsA * misc::asinh(tmp2);
  return result;
}

template <typename Config>
void RateAndStateFastVelocityInitializer<Config>::addAdditionalParameters(
    std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateInitializer<Config>::addAdditionalParameters(parameterToStorageMap, dynRup, it);
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSRateAndStateFastVelocityWeakening const* const>(
          dynRup);
  typename Config::RealT(*rsSrW)[misc::numPaddedPoints<Config>] = it->var(concreteLts->rsSrW);
  parameterToStorageMap.insert({"rs_srW", (typename Config::RealT*)rsSrW});
}

template <typename Config>
void RateAndStateThermalPressurizationInitializer<Config>::initializeFault(
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  RateAndStateInitializer<Config>::initializeFault(dynRup, dynRupTree);

  auto* concreteLts = dynamic_cast<
      seissol::initializers::LTSRateAndStateThermalPressurization<Config> const* const>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    typename Config::RealT(*temperature)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->temperature);
    typename Config::RealT(*pressure)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->pressure);
    typename Config::RealT(*theta)[misc::numPaddedPoints<Config>][misc::numberOfTPGridPoints] =
        it->var(concreteLts->theta);
    typename Config::RealT(*sigma)[misc::numPaddedPoints<Config>][misc::numberOfTPGridPoints] =
        it->var(concreteLts->sigma);
    typename Config::RealT(
        *thetaTmpBuffer)[misc::numPaddedPoints<Config>][misc::numberOfTPGridPoints] =
        it->var(concreteLts->thetaTmpBuffer);
    typename Config::RealT(
        *sigmaTmpBuffer)[misc::numPaddedPoints<Config>][misc::numberOfTPGridPoints] =
        it->var(concreteLts->sigmaTmpBuffer);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; ++pointIndex) {
        temperature[ltsFace][pointIndex] = drParameters->initialTemperature;
        pressure[ltsFace][pointIndex] = drParameters->initialPressure;
        for (unsigned tpGridPointIndex = 0; tpGridPointIndex < misc::numberOfTPGridPoints;
             ++tpGridPointIndex) {
          theta[ltsFace][pointIndex][tpGridPointIndex] = 0.0;
          sigma[ltsFace][pointIndex][tpGridPointIndex] = 0.0;
          thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = 0.0;
          sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = 0.0;
        }
      }
    }
  }
}

template <typename Config>
void RateAndStateThermalPressurizationInitializer<Config>::addAdditionalParameters(
    std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateFastVelocityInitializer<Config>::addAdditionalParameters(
      parameterToStorageMap, dynRup, it);

  auto* concreteLts = dynamic_cast<
      seissol::initializers::LTSRateAndStateThermalPressurization<Config> const* const>(dynRup);

  typename Config::RealT(*halfWidthShearZone)[misc::numPaddedPoints<Config>] =
      it->var(concreteLts->halfWidthShearZone);
  typename Config::RealT(*hydraulicDiffusivity)[misc::numPaddedPoints<Config>] =
      it->var(concreteLts->hydraulicDiffusivity);
  parameterToStorageMap.insert(
      {"tp_halfWidthShearZone", (typename Config::RealT*)halfWidthShearZone});
  parameterToStorageMap.insert(
      {"tp_hydraulicDiffusivity", (typename Config::RealT*)hydraulicDiffusivity});
}
} // namespace seissol::dr::initializers
