#include "DynamicRupture/Output/OutputAux.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Numerical_aux/BasisFunction.h"
#include "ReceiverBasedOutput.hpp"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"
#include <unordered_map>

namespace seissol::dr::output {
void ReceiverBasedOutput::setLtsData(seissol::initializers::LTSTree* userWpTree,
                                     seissol::initializers::LTS* userWpDescr,
                                     seissol::initializers::Lut* userWpLut,
                                     seissol::initializers::LTSTree* userDrTree,
                                     seissol::initializers::DynamicRupture* userDrDescr) {
  wpTree = userWpTree;
  wpDescr = userWpDescr;
  wpLut = userWpLut;
  drTree = userDrTree;
  drDescr = userDrDescr;
}

void ReceiverBasedOutput::getDofs(real dofs[tensor::Q::size()], int meshId) {
  // get DOFs from 0th derivatives
  assert((wpLut->lookup(wpDescr->cellInformation, meshId).ltsSetup >> 9) % 2 == 1);

  real* derivatives = wpLut->lookup(wpDescr->derivatives, meshId);
  std::copy(&derivatives[0], &derivatives[tensor::dQ::Size[0]], &dofs[0]);
}

void ReceiverBasedOutput::getNeighbourDofs(real dofs[tensor::Q::size()], int meshId, int side) {
  real* derivatives = wpLut->lookup(wpDescr->faceNeighbors, meshId)[side];
  assert(derivatives != nullptr);

  std::copy(&derivatives[0], &derivatives[tensor::dQ::Size[0]], &dofs[0]);
}

void ReceiverBasedOutput::calcFaultOutput(const OutputType type,
                                          ReceiverBasedOutputData& outputData,
                                          const GeneralParamsT& generalParams,
                                          double time) {

  const size_t level = (type == OutputType::AtPickpoint) ? outputData.currentCacheLevel : 0;
  const auto faultInfos = meshReader->getFault();

  for (size_t i = 0; i < outputData.receiverPoints.size(); ++i) {

    assert(outputData.receiverPoints[i].isInside == true &&
           "A receiver must be inside of a fault. Error in pre-processing");

    const auto faceIndex = outputData.receiverPoints[i].faultFaceIndex;
    assert(faceIndex != -1 && "receiver is not initialized");
    local = LocalInfo{};

    const auto ltsMap = (*faceToLtsMap)[faceIndex];
    local.layer = ltsMap.first;
    local.ltsId = ltsMap.second;
    local.nearestGpIndex = outputData.receiverPoints[i].nearestGpIndex;
    local.nearestInternalGpIndex = outputData.receiverPoints[i].nearestInternalGpIndex;

    local.waveSpeedsPlus = &((local.layer->var(drDescr->waveSpeedsPlus))[local.ltsId]);
    local.waveSpeedsMinus = &((local.layer->var(drDescr->waveSpeedsMinus))[local.ltsId]);

    const auto faultInfo = faultInfos[faceIndex];

    real dofsPlus[tensor::Q::size()]{};
    getDofs(dofsPlus, faultInfo.element);

    real dofsMinus[tensor::Q::size()]{};
    if (faultInfo.neighborElement >= 0) {
      getDofs(dofsMinus, faultInfo.neighborElement);
    } else {
      getNeighbourDofs(dofsMinus, faultInfo.element, faultInfo.side);
    }

    const auto* initStresses = local.layer->var(drDescr->initialStressInFaultCS);
    const auto* initStress = initStresses[local.ltsId][local.nearestGpIndex];

    local.frictionCoefficient = (local.layer->var(drDescr->mu))[local.ltsId][local.nearestGpIndex];
    local.stateVariable = this->computeStateVariable();

    local.iniTraction1 = initStress[3];
    local.iniTraction2 = initStress[5];
    local.iniNormalTraction = initStress[0];
    local.fluidPressure = this->computeFluidPressure();

    const auto* const normal = outputData.faultDirections[i].faceNormal;
    const auto* const tangent1 = outputData.faultDirections[i].tangent1;
    const auto* const tangent2 = outputData.faultDirections[i].tangent2;
    const auto* strike = outputData.faultDirections[i].strike;
    const auto* dip = outputData.faultDirections[i].dip;

    auto* phiPlusSide = outputData.basisFunctions[i].plusSide.data();
    auto* phiMinusSide = outputData.basisFunctions[i].minusSide.data();

    seissol::dynamicRupture::kernel::evaluateFaceAlignedDOFSAtPoint kernel;
    kernel.Tinv = outputData.glbToFaceAlignedData[i].data();

    kernel.Q = dofsPlus;
    kernel.basisFunctionsAtPoint = phiPlusSide;
    kernel.QAtPoint = local.faceAlignedValuesPlus;
    kernel.execute();

    kernel.Q = dofsMinus;
    kernel.basisFunctionsAtPoint = phiMinusSide;
    kernel.QAtPoint = local.faceAlignedValuesMinus;
    kernel.execute();

    this->computeLocalStresses();
    const real strength = this->computeLocalStrength();
    this->updateLocalTractions(strength);

    seissol::dynamicRupture::kernel::rotateInitStress alignAlongDipAndStrikeKernel;
    alignAlongDipAndStrikeKernel.stressRotationMatrix =
        outputData.stressGlbToDipStrikeAligned[i].data();
    alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
        outputData.stressFaceAlignedToGlb[i].data();

    std::array<real, 6> tmpVector{};
    tmpVector[0] = local.transientNormalTraction;
    tmpVector[1] = local.faceAlignedStress22;
    tmpVector[2] = local.faceAlignedStress33;
    tmpVector[3] = local.updatedTraction1;
    tmpVector[4] = local.faceAlignedStress23;
    tmpVector[5] = local.updatedTraction2;

    alignAlongDipAndStrikeKernel.initialStress = tmpVector.data();
    std::array<real, 6> rotatedUpdatedStress{};
    alignAlongDipAndStrikeKernel.rotatedStress = rotatedUpdatedStress.data();
    alignAlongDipAndStrikeKernel.execute();

    tmpVector[0] = local.transientNormalTraction;
    tmpVector[1] = local.faceAlignedStress22;
    tmpVector[2] = local.faceAlignedStress33;
    tmpVector[3] = local.faceAlignedStress12;
    tmpVector[4] = local.faceAlignedStress23;
    tmpVector[5] = local.faceAlignedStress13;

    alignAlongDipAndStrikeKernel.initialStress = tmpVector.data();
    std::array<real, 6> rotatedStress{};
    alignAlongDipAndStrikeKernel.rotatedStress = rotatedStress.data();
    alignAlongDipAndStrikeKernel.execute();

    switch (generalParams.slipRateOutputType) {
    case SlipRateOutputType::TractionsAndFailure: {
      this->computeSlipRate(rotatedUpdatedStress, rotatedStress);
      break;
    }
    case SlipRateOutputType::VelocityDifference: {
      this->computeSlipRate(tangent1, tangent2, strike, dip);
      break;
    }
    }

    adjustRotatedUpdatedStress(rotatedUpdatedStress, rotatedStress);

    auto& slipRate = std::get<VariableID::SlipRate>(outputData.vars);
    if (slipRate.isActive) {
      slipRate(DirectionID::Strike, level, i) = local.slipRateStrike;
      slipRate(DirectionID::Dip, level, i) = local.slipRateDip;
    }

    auto& transientTractions = std::get<VariableID::TransientTractions>(outputData.vars);
    if (transientTractions.isActive) {
      transientTractions(DirectionID::Strike, level, i) = rotatedUpdatedStress[3];
      transientTractions(DirectionID::Dip, level, i) = rotatedUpdatedStress[5];
      transientTractions(DirectionID::Normal, level, i) =
          local.transientNormalTraction - local.fluidPressure;
    }

    auto& frictionAndState = std::get<VariableID::FrictionAndState>(outputData.vars);
    if (frictionAndState.isActive) {
      frictionAndState(ParamID::FrictionCoefficient, level, i) = local.frictionCoefficient;
      frictionAndState(ParamID::State, level, i) = local.stateVariable;
    }

    auto& ruptureTime = std::get<VariableID::RuptureTime>(outputData.vars);
    if (ruptureTime.isActive) {
      auto* rt = local.layer->var(drDescr->ruptureTime);
      ruptureTime(level, i) = rt[local.ltsId][local.nearestGpIndex];
    }

    auto& normalVelocity = std::get<VariableID::NormalVelocity>(outputData.vars);
    if (normalVelocity.isActive) {
      normalVelocity(level, i) = local.faultNormalVelocity;
    }

    auto& accumulatedSlip = std::get<VariableID::AccumulatedSlip>(outputData.vars);
    if (accumulatedSlip.isActive) {
      auto* slip = local.layer->var(drDescr->accumulatedSlipMagnitude);
      accumulatedSlip(level, i) = slip[local.ltsId][local.nearestGpIndex];
    }

    auto& totalTractions = std::get<VariableID::TotalTractions>(outputData.vars);
    if (totalTractions.isActive) {
      std::array<real, tensor::rotatedStress::size()> rotatedInitStress{};
      alignAlongDipAndStrikeKernel.initialStress = initStress;
      alignAlongDipAndStrikeKernel.rotatedStress = rotatedInitStress.data();
      alignAlongDipAndStrikeKernel.execute();

      totalTractions(DirectionID::Strike, level, i) =
          rotatedUpdatedStress[3] + rotatedInitStress[3];
      totalTractions(DirectionID::Dip, level, i) = rotatedUpdatedStress[5] + rotatedInitStress[5];
      totalTractions(DirectionID::Normal, level, i) =
          local.transientNormalTraction - local.fluidPressure + rotatedInitStress[0];
    }

    auto& ruptureVelocity = std::get<VariableID::RuptureVelocity>(outputData.vars);
    if (ruptureVelocity.isActive) {
      auto& jacobiT2d = outputData.jacobianT2d[i];
      ruptureVelocity(level, i) = this->computeRuptureVelocity(jacobiT2d);
    }

    auto& peakSlipsRate = std::get<VariableID::PeakSlipRate>(outputData.vars);
    if (peakSlipsRate.isActive) {
      auto* peakSR = local.layer->var(drDescr->peakSlipRate);
      peakSlipsRate(level, i) = peakSR[local.ltsId][local.nearestGpIndex];
    }

    auto& dynamicStressTime = std::get<VariableID::DynamicStressTime>(outputData.vars);
    if (dynamicStressTime.isActive) {
      auto* dynStressTime = (local.layer->var(drDescr->dynStressTime));
      dynamicStressTime(level, i) = dynStressTime[local.ltsId][local.nearestGpIndex];
    }

    auto& slipVectors = std::get<VariableID::Slip>(outputData.vars);
    if (slipVectors.isActive) {
      VrtxCoords crossProduct = {0.0, 0.0, 0.0};
      MeshTools::cross(strike, tangent1, crossProduct);

      const double cos1 = MeshTools::dot(strike, tangent1);
      const double scalarProd = MeshTools::dot(crossProduct, normal);

      // Note: cos1**2 can be greater than 1.0 because of rounding errors -> min
      double sin1 = std::sqrt(1.0 - std::min(1.0, cos1 * cos1));
      sin1 = (scalarProd > 0) ? sin1 : -sin1;

      auto* slip1 = local.layer->var(drDescr->slip1);
      auto* slip2 = local.layer->var(drDescr->slip2);

      slipVectors(DirectionID::Strike, level, i) = cos1 * slip1[local.ltsId][local.nearestGpIndex] -
                                                   sin1 * slip2[local.ltsId][local.nearestGpIndex];

      slipVectors(DirectionID::Dip, level, i) = sin1 * slip1[local.ltsId][local.nearestGpIndex] +
                                                cos1 * slip2[local.ltsId][local.nearestGpIndex];
    }
    this->outputSpecifics(outputData, level, i);
  }

  if (type == OutputType::AtPickpoint) {
    outputData.cachedTime[outputData.currentCacheLevel] = time;
    outputData.currentCacheLevel += 1;
  }
}

void ReceiverBasedOutput::computeLocalStresses() {
  const auto& impAndEta = ((local.layer->var(drDescr->impAndEta))[local.ltsId]);
  const real normalDivisor = 1.0 / (impAndEta.zpNeig + impAndEta.zp);
  const real shearDivisor = 1.0 / (impAndEta.zsNeig + impAndEta.zs);

  auto diff = [this](int i) {
    return this->local.faceAlignedValuesMinus[i] - this->local.faceAlignedValuesPlus[i];
  };

  local.faceAlignedStress12 =
      local.faceAlignedValuesPlus[3] +
      ((diff(3) + impAndEta.zsNeig * diff(7)) * impAndEta.zs) * shearDivisor;

  local.faceAlignedStress13 =
      local.faceAlignedValuesPlus[5] +
      ((diff(5) + impAndEta.zsNeig * diff(8)) * impAndEta.zs) * shearDivisor;

  local.transientNormalTraction =
      local.faceAlignedValuesPlus[0] +
      ((diff(0) + impAndEta.zpNeig * diff(6)) * impAndEta.zp) * normalDivisor;

  local.faultNormalVelocity =
      local.faceAlignedValuesPlus[6] +
      (local.transientNormalTraction - local.faceAlignedValuesPlus[0]) * impAndEta.invZp;

  real missingSigmaValues = (local.transientNormalTraction - local.faceAlignedValuesPlus[0]);
  missingSigmaValues *= (1.0 - 2.0 * std::pow(local.waveSpeedsPlus->sWaveVelocity /
                                                  local.waveSpeedsPlus->pWaveVelocity,
                                              2));

  local.faceAlignedStress22 = local.faceAlignedValuesPlus[1] + missingSigmaValues;
  local.faceAlignedStress33 = local.faceAlignedValuesPlus[2] + missingSigmaValues;
  local.faceAlignedStress23 = local.faceAlignedValuesPlus[4];
}

void ReceiverBasedOutput::updateLocalTractions(real strength) {
  const auto component1 = local.iniTraction1 + local.faceAlignedStress12;
  const auto component2 = local.iniTraction2 + local.faceAlignedStress13;
  const auto tracEla = misc::magnitude(component1, component2);

  if (tracEla > std::abs(strength)) {
    local.updatedTraction1 =
        ((local.iniTraction1 + local.faceAlignedStress12) / tracEla) * strength;
    local.updatedTraction2 =
        ((local.iniTraction2 + local.faceAlignedStress13) / tracEla) * strength;

    // update stress change
    local.updatedTraction1 -= local.iniTraction1;
    local.updatedTraction2 -= local.iniTraction2;
  } else {
    local.updatedTraction1 = local.faceAlignedStress12;
    local.updatedTraction2 = local.faceAlignedStress13;
  }
}

void ReceiverBasedOutput::computeSlipRate(std::array<real, 6>& rotatedUpdatedStress,
                                          std::array<real, 6>& rotatedStress) {

  const auto& impAndEta = ((local.layer->var(drDescr->impAndEta))[local.ltsId]);
  local.slipRateStrike = -impAndEta.invEtaS * (rotatedUpdatedStress[3] - rotatedStress[3]);
  local.slipRateDip = -impAndEta.invEtaS * (rotatedUpdatedStress[5] - rotatedStress[5]);
}

void ReceiverBasedOutput::computeSlipRate(const double* tangent1,
                                          const double* tangent2,
                                          const double* strike,
                                          const double* dip) {
  local.slipRateStrike = static_cast<real>(0.0);
  local.slipRateDip = static_cast<real>(0.0);

  for (size_t i = 0; i < 3; ++i) {
    const real factorMinus = (local.faceAlignedValuesMinus[7] * tangent1[i] +
                              local.faceAlignedValuesMinus[8] * tangent2[i]);

    const real factorPlus = (local.faceAlignedValuesPlus[7] * tangent1[i] +
                             local.faceAlignedValuesPlus[8] * tangent2[i]);

    local.slipRateStrike += (factorMinus - factorPlus) * strike[i];
    local.slipRateDip += (factorMinus - factorPlus) * dip[i];
  }
}

real ReceiverBasedOutput::computeRuptureVelocity(Eigen::Matrix<real, 2, 2>& jacobiT2d) {
  auto* ruptureTime = (local.layer->var(drDescr->ruptureTime))[local.ltsId];
  real ruptureVelocity = 0.0;

  bool needsUpdate{true};
  for (size_t point = 0; point < misc::numberOfBoundaryGaussPoints; ++point) {
    if (ruptureTime[point] == 0.0) {
      needsUpdate = false;
    }
  }

  if (needsUpdate) {
    constexpr int numPoly = CONVERGENCE_ORDER - 1;
    constexpr int numDegFr2d = (numPoly + 1) * (numPoly + 2) / 2;
    std::array<double, numDegFr2d> projectedRT{};
    projectedRT.fill(0.0);

    std::array<double, 2 * numDegFr2d> phiAtPoint{};
    phiAtPoint.fill(0.0);

    auto chiTau2dPoints =
        init::quadpoints::view::create(const_cast<real*>(init::quadpoints::Values));
    auto weights = init::quadweights::view::create(const_cast<real*>(init::quadweights::Values));

    auto* rt = local.layer->var(drDescr->ruptureTime);
    for (size_t jBndGP = 0; jBndGP < misc::numberOfBoundaryGaussPoints; ++jBndGP) {
      real chi = chiTau2dPoints(jBndGP, 0);
      real tau = chiTau2dPoints(jBndGP, 1);
      basisFunction::TriDubiner::evaluatePolynomials(phiAtPoint.data(), chi, tau, numPoly);

      for (size_t d = 0; d < numDegFr2d; ++d) {
        projectedRT[d] += weights(jBndGP) * rt[local.ltsId][jBndGP] * phiAtPoint[d];
      }
    }

    auto m2inv =
        seissol::init::M2inv::view::create(const_cast<real*>(seissol::init::M2inv::Values));
    for (size_t d = 0; d < numDegFr2d; ++d) {
      projectedRT[d] *= m2inv(d, d);
    }

    const real chi = chiTau2dPoints(local.nearestInternalGpIndex, 0);
    const real tau = chiTau2dPoints(local.nearestInternalGpIndex, 1);

    basisFunction::TriDubiner::evaluateGradPolynomials(phiAtPoint.data(), chi, tau, numPoly);

    real dTdChi{0.0};
    real dTdTau{0.0};
    for (size_t d = 0; d < numDegFr2d; ++d) {
      dTdChi += projectedRT[d] * phiAtPoint[2 * d];
      dTdTau += projectedRT[d] * phiAtPoint[2 * d + 1];
    }
    const real dTdX = jacobiT2d(0, 0) * dTdChi + jacobiT2d(0, 1) * dTdTau;
    const real dTdY = jacobiT2d(1, 0) * dTdChi + jacobiT2d(1, 1) * dTdTau;

    const real slowness = misc::magnitude(dTdX, dTdY);
    ruptureVelocity = (slowness == 0.0) ? 0.0 : 1.0 / slowness;
  }

  return ruptureVelocity;
}
} // namespace seissol::dr::output
