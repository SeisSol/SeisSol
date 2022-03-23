#include "ImposedSlipRatesInitializer.h"

#include <utils/logger.h>

namespace seissol::dr::initializers {
void ImposedSlipRatesInitializer::initializeFault(seissol::initializers::DynamicRupture* dynRup,
                                                  seissol::initializers::LTSTree* dynRupTree,
                                                  seissol::Interoperability* eInteroperability) {
  logInfo() << "Initializing Fault, using a quadrature rule with "
            << misc::numberOfBoundaryGaussPoints << " points.";
  seissol::initializers::FaultParameterDB faultParameterDB;

  for (auto it = dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, real*> parameterToStorageMap;

    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRates*>(dynRup);
    real(*strikeSlip)[misc::numPaddedPoints] = it->var(concreteLts->strikeSlip);
    real(*dipSlip)[misc::numPaddedPoints] = it->var(concreteLts->dipSlip);
    real(*onsetTime)[misc::numPaddedPoints] = it->var(concreteLts->onsetTime);
    real(*tauS)[misc::numPaddedPoints] = it->var(concreteLts->tauS);
    real(*tauR)[misc::numPaddedPoints] = it->var(concreteLts->tauR);
    parameterToStorageMap.insert({"strike_slip", (real*)strikeSlip});
    parameterToStorageMap.insert({"dip_slip", (real*)dipSlip});
    parameterToStorageMap.insert({"rupture_onset", (real*)onsetTime});
    parameterToStorageMap.insert({"tau_S", (real*)tauS});
    parameterToStorageMap.insert({"tau_R", (real*)tauR});

    // read parameters from yaml file
    for (const auto& parameterStoragePair : parameterToStorageMap) {
      faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
    }
    const auto faceIDs = getFaceIDsInIterator(dynRup, it);
    queryModel(faultParameterDB, faceIDs);

    real(*nucleationStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(dynRup->nucleationStressInFaultCS);
    real(*initialStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(dynRup->initialStressInFaultCS);

    // Set initial and nucleation stress to zero, these are not needed for this FL
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        for (unsigned int dim = 0; dim < 6; ++dim) {
          initialStressInFaultCS[ltsFace][pointIndex][dim] = 0;
          nucleationStressInFaultCS[ltsFace][pointIndex][dim] = 0;
        }
      }
    }

    // ensure that tauR is larger than tauS and that tauS and tauR are greater than 0 (the contrary
    // can happen due to ASAGI interpolation)
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        tauS[ltsFace][pointIndex] = std::max(static_cast<real>(0.0), tauS[ltsFace][pointIndex]);
        tauR[ltsFace][pointIndex] = std::max(tauR[ltsFace][pointIndex], tauS[ltsFace][pointIndex]);
      }
    }
  }
}
} // namespace seissol::dr::initializers