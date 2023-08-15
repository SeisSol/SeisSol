#include "ImposedSlipRatesInitializer.h"

#include "Model/common.hpp"
#include "SeisSol.h"
#include <utils/logger.h>

namespace seissol::dr::initializers {
template <typename Config>
void ImposedSlipRatesInitializer<Config>::initializeFault(
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  const int rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Initializing Fault, using a quadrature rule with "
                << misc::numberOfBoundaryGaussPoints<Config> << " points.";
  seissol::initializers::FaultParameterDB faultParameterDB;

  for (auto it = dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    // parameters to be read from fault parameters yaml file
    std::unordered_map<std::string, typename Config::RealT*> parameterToStorageMap;

    auto* concreteLts =
        dynamic_cast<seissol::initializers::LTSImposedSlipRates<Config> const* const>(dynRup);
    auto* imposedSlipDirection1 = it->var(concreteLts->imposedSlipDirection1);
    auto* imposedSlipDirection2 = it->var(concreteLts->imposedSlipDirection2);
    auto* onsetTime = it->var(concreteLts->onsetTime);

    // First read slip in strike/dip direction. Later we will rotate this to the face aligned
    // coordinate system.
    using VectorOfArraysT =
        std::vector<std::array<typename Config::RealT, misc::numPaddedPoints<Config>>>;
    VectorOfArraysT strikeSlip(it->getNumberOfCells());
    VectorOfArraysT dipSlip(it->getNumberOfCells());
    parameterToStorageMap.insert({"strike_slip", strikeSlip.data()->data()});
    parameterToStorageMap.insert({"dip_slip", dipSlip.data()->data()});
    parameterToStorageMap.insert({"rupture_onset", (typename Config::RealT*)onsetTime});

    // get additional parameters (for derived friction laws)
    addAdditionalParameters(parameterToStorageMap, dynRup, it);

    // read parameters from yaml file
    for (const auto& parameterStoragePair : parameterToStorageMap) {
      faultParameterDB.addParameter(parameterStoragePair.first, parameterStoragePair.second);
    }
    const auto faceIDs = getFaceIDsInIterator(dynRup, it);
    queryModel(faultParameterDB, faceIDs);

    rotateSlipToFaultCS(
        dynRup, it, strikeSlip, dipSlip, imposedSlipDirection1, imposedSlipDirection2);

    typename Config::RealT(*nucleationStressInFaultCS)[misc::numPaddedPoints<Config>][6] =
        it->var(dynRup->nucleationStressInFaultCS);
    typename Config::RealT(*initialStressInFaultCS)[misc::numPaddedPoints<Config>][6] =
        it->var(dynRup->initialStressInFaultCS);

    // Set initial and nucleation stress to zero, these are not needed for this FL
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; ++pointIndex) {
        for (unsigned int dim = 0; dim < 6; ++dim) {
          initialStressInFaultCS[ltsFace][pointIndex][dim] = 0;
          nucleationStressInFaultCS[ltsFace][pointIndex][dim] = 0;
        }
      }
    }

    fixInterpolatedSTFParameters(dynRup, it);

    initializeOtherVariables(dynRup, it);
  }
}

template <typename Config>
void ImposedSlipRatesInitializer<Config>::rotateSlipToFaultCS(
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSTree::leaf_iterator& it,
    std::vector<std::array<typename Config::RealT, misc::numPaddedPoints<Config>>> const&
        strikeSlip,
    std::vector<std::array<typename Config::RealT, misc::numPaddedPoints<Config>>> const& dipSlip,
    typename Config::RealT (*imposedSlipDirection1)[misc::numPaddedPoints<Config>],
    typename Config::RealT (*imposedSlipDirection2)[misc::numPaddedPoints<Config>]) {
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    const auto& drFaceInformation = it->var(dynRup->faceInformation);
    const unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
    const Fault& fault = seissol::SeisSol::main.meshReader().getFault().at(meshFace);

    VrtxCoords strikeVector{};
    VrtxCoords dipVector{};
    misc::computeStrikeAndDipVectors(fault.normal, strikeVector, dipVector);

    // cos^2 can be greater than 1 because of rounding errors
    typename Config::RealT cos =
        std::clamp(MeshTools::dot(strikeVector, fault.tangent1), -1.0, 1.0);
    VrtxCoords crossProduct{};
    MeshTools::cross(strikeVector, fault.tangent1, crossProduct);
    typename Config::RealT scalarProduct = MeshTools::dot(crossProduct, fault.normal);
    typename Config::RealT sin = std::sqrt(1 - cos * cos) * std::copysign(1.0, scalarProduct);
    for (size_t pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; ++pointIndex) {
      imposedSlipDirection1[ltsFace][pointIndex] =
          cos * strikeSlip[ltsFace][pointIndex] + sin * dipSlip[ltsFace][pointIndex];
      imposedSlipDirection2[ltsFace][pointIndex] =
          -sin * strikeSlip[ltsFace][pointIndex] + cos * dipSlip[ltsFace][pointIndex];
    }
  }
}

template <typename Config>
void ImposedSlipRatesInitializer<Config>::fixInterpolatedSTFParameters(
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  // do nothing
}

template <typename Config>
void ImposedSlipRatesYoffeInitializer<Config>::addAdditionalParameters(
    std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSImposedSlipRatesYoffe<Config> const* const>(dynRup);
  typename Config::RealT(*tauS)[misc::numPaddedPoints<Config>] = it->var(concreteLts->tauS);
  typename Config::RealT(*tauR)[misc::numPaddedPoints<Config>] = it->var(concreteLts->tauR);
  parameterToStorageMap.insert({"tau_S", (typename Config::RealT*)tauS});
  parameterToStorageMap.insert({"tau_R", (typename Config::RealT*)tauR});
}

template <typename Config>
void ImposedSlipRatesYoffeInitializer<Config>::fixInterpolatedSTFParameters(
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSImposedSlipRatesYoffe<Config> const* const>(dynRup);
  typename Config::RealT(*tauS)[misc::numPaddedPoints<Config>] = it->var(concreteLts->tauS);
  typename Config::RealT(*tauR)[misc::numPaddedPoints<Config>] = it->var(concreteLts->tauR);
  // ensure that tauR is larger than tauS and that tauS and tauR are greater than 0 (the contrary
  // can happen due to ASAGI interpolation)
  for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
    for (unsigned int pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; ++pointIndex) {
      tauS[ltsFace][pointIndex] =
          std::max(static_cast<typename Config::RealT>(0.0), tauS[ltsFace][pointIndex]);
      tauR[ltsFace][pointIndex] = std::max(tauR[ltsFace][pointIndex], tauS[ltsFace][pointIndex]);
    }
  }
}

template <typename Config>
void ImposedSlipRatesGaussianInitializer<Config>::addAdditionalParameters(
    std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSImposedSlipRatesGaussian<Config> const* const>(dynRup);
  typename Config::RealT(*riseTime)[misc::numPaddedPoints<Config>] = it->var(concreteLts->riseTime);
  parameterToStorageMap.insert({"rupture_rise_time", (typename Config::RealT*)riseTime});
}
} // namespace seissol::dr::initializers
