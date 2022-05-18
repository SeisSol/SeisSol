#include "ReceiverBasedOutput.hpp"
#include "DynamicRupture/Output/OutputAux.hpp"
#include "Initializer/tree/Layer.hpp"
#include "generated_code/tensor.h"
#include "generated_code/kernel.h"
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
  auto ltsSetup = wpLut->lookup(wpDescr->cellInformation, meshId).ltsSetup;
  assert((ltsSetup >> 9) % 2 == 1);

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

  size_t level = (type == OutputType::AtPickpoint) ? outputData.currentCacheLevel : 0;
  auto faultInfos = meshReader->getFault();
  for (size_t i = 0; i < outputData.receiverPoints.size(); ++i) {

    assert(outputData.receiverPoints[i].isInside == true &&
           "A receiver must be inside of a fault. Error in pre-processing");

    const auto faceIndex = outputData.receiverPoints[i].faultFaceIndex;
    assert(faceIndex != -1 && "receiver is not initialized");
    local = LocalInfo{};

    auto ltsMap = (*faceToLtsMap)[faceIndex];
    local.layer = ltsMap.first;
    local.ltsId = ltsMap.second;
    local.nearestGpIndex = outputData.receiverPoints[i].nearestGpIndex;

    local.waveSpeedsPlus = &((local.layer->var(drDescr->waveSpeedsPlus))[local.ltsId]);
    local.waveSpeedsMinus = &((local.layer->var(drDescr->waveSpeedsMinus))[local.ltsId]);

    auto faultInfo = faultInfos[faceIndex];

    real dofsPlus[tensor::Q::size()]{};
    getDofs(dofsPlus, faultInfo.element);

    real dofsMinus[tensor::Q::size()]{};
    if (faultInfo.neighborElement >= 0) {
      getDofs(dofsMinus, faultInfo.neighborElement);
    } else {
      getNeighbourDofs(dofsMinus, faultInfo.element, faultInfo.side);
    }

    auto* initStresses = local.layer->var(drDescr->initialStressInFaultCS);
    auto* initStress = initStresses[local.ltsId][local.nearestGpIndex];

    local.mu = (local.layer->var(drDescr->mu))[local.ltsId][local.nearestGpIndex];
    local.sXY = initStress[3];
    local.sXZ = initStress[5];
    local.p0 = initStress[0];
    local.pf = this->computePf();

    const auto* const normal = outputData.faultDirections[i].faceNormal;
    const auto* const tangent1 = outputData.faultDirections[i].tangent1;
    const auto* const tangent2 = outputData.faultDirections[i].tangent2;
    auto* const strike = outputData.faultDirections[i].strike;
    const auto* dip = outputData.faultDirections[i].dip;

    auto* phiPlusSide = outputData.basisFunctions[i].plusSide.data();
    auto* phiMinusSide = outputData.basisFunctions[i].minusSide.data();

    seissol::dynamicRupture::kernel::computeFaceAlignedValues kernel;
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
    real strength = this->computeLocalStrength();
    this->computeLocalTraction(strength);

    seissol::dynamicRupture::kernel::rotateInitStress alignAlongDipAndStrikeKernel;
    alignAlongDipAndStrikeKernel.stressRotationMatrix =
        outputData.stressGlbToDipStrikeAligned[i].data();
    alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
        outputData.stressFaceAlignedToGlb[i].data();

    std::array<real, 6> tmpVector{};
    tmpVector[0] = local.p;
    tmpVector[1] = local.yyStress;
    tmpVector[2] = local.zzStress;
    tmpVector[3] = local.xyTraction;
    tmpVector[4] = local.yzStress;
    tmpVector[5] = local.xzTraction;

    alignAlongDipAndStrikeKernel.initialStress = tmpVector.data();
    std::array<real, 6> rotatedTraction{};
    alignAlongDipAndStrikeKernel.rotatedStress = rotatedTraction.data();
    alignAlongDipAndStrikeKernel.execute();

    tmpVector[0] = local.p;
    tmpVector[1] = local.yyStress;
    tmpVector[2] = local.zzStress;
    tmpVector[3] = local.xyStress;
    tmpVector[4] = local.yzStress;
    tmpVector[5] = local.xzStress;

    alignAlongDipAndStrikeKernel.initialStress = tmpVector.data();
    std::array<real, 6> rotatedLocalStress{};
    alignAlongDipAndStrikeKernel.rotatedStress = rotatedLocalStress.data();
    alignAlongDipAndStrikeKernel.execute();

    if (generalParams.slipRateOutputType) {
      this->computeSlipAndRate(rotatedTraction, rotatedLocalStress);
    } else {
      this->computeSlipAndRate(tangent1, tangent2, strike, dip);
    }

    adjustRotatedTractionAndStresses(rotatedTraction, rotatedLocalStress);

    auto& slipRate = std::get<VariableID::SlipRate>(outputData.vars);
    if (slipRate.isActive) {
      slipRate(DirectionID::Strike, level, i) = local.srS;
      slipRate(DirectionID::Dip, level, i) = local.srD;
    }

    auto& transientTractions = std::get<VariableID::TransientTractions>(outputData.vars);
    if (transientTractions.isActive) {
      transientTractions(DirectionID::Strike, level, i) = rotatedTraction[3];
      transientTractions(DirectionID::Dip, level, i) = rotatedTraction[5];
      transientTractions(DirectionID::Normal, level, i) = local.p - local.pf;
    }

    auto& frictionAndState = std::get<VariableID::FrictionAndState>(outputData.vars);
    if (frictionAndState.isActive) {
      frictionAndState(ParamID::FrictionCoefficient, level, i) = local.mu;
    }

    auto& ruptureTime = std::get<VariableID::RuptureTime>(outputData.vars);
    if (ruptureTime.isActive) {
      auto* rt = local.layer->var(drDescr->ruptureTime);
      ruptureTime(level, i) = rt[local.ltsId][local.nearestGpIndex];
    }

    auto& normalVelocity = std::get<VariableID::NormalVelocity>(outputData.vars);
    if (normalVelocity.isActive) {
      normalVelocity(level, i) = local.u;
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

      totalTractions(DirectionID::Strike, level, i) = rotatedTraction[3] + rotatedInitStress[3];
      totalTractions(DirectionID::Dip, level, i) = rotatedTraction[5] + rotatedInitStress[5];
      totalTractions(DirectionID::Normal, level, i) = local.p - local.pf + rotatedInitStress[0];
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

      double cos1 = MeshTools::dot(strike, tangent1);
      double scalarProd = MeshTools::dot(crossProduct, normal);

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
  auto& impAndEta = ((local.layer->var(drDescr->impAndEta))[local.ltsId]);
  real norDivisor = 1.0 / (impAndEta.zpNeig + impAndEta.zp);
  real shearDivisor = 1.0 / (impAndEta.zsNeig + impAndEta.zs);

  auto diff = [this](int i) {
    return this->local.faceAlignedValuesMinus[i] - this->local.faceAlignedValuesPlus[i];
  };

  local.xyStress = local.faceAlignedValuesPlus[3] +
                   ((diff(3) + impAndEta.zsNeig * diff(7)) * impAndEta.zs) * shearDivisor;

  local.xzStress = local.faceAlignedValuesPlus[5] +
                   ((diff(5) + impAndEta.zsNeig * diff(8)) * impAndEta.zs) * shearDivisor;

  local.p = local.faceAlignedValuesPlus[0] +
            ((diff(0) + impAndEta.zpNeig * diff(6)) * impAndEta.zp) * norDivisor;

  local.u =
      local.faceAlignedValuesPlus[6] + (local.p - local.faceAlignedValuesPlus[0]) * impAndEta.invZp;

  real missingSigmaValues = (local.p - local.faceAlignedValuesPlus[0]);
  missingSigmaValues *= (1.0 - 2.0 * std::pow(local.waveSpeedsPlus->sWaveVelocity /
                                                  local.waveSpeedsPlus->pWaveVelocity,
                                              2));

  local.yyStress = local.faceAlignedValuesPlus[1] + missingSigmaValues;
  local.zzStress = local.faceAlignedValuesPlus[2] + missingSigmaValues;
  local.yzStress = local.faceAlignedValuesPlus[4];

  local.tracEla =
      std::sqrt(std::pow(local.sXY + local.xyStress, 2) + std::pow(local.sXZ + local.xzStress, 2));
}

void ReceiverBasedOutput::computeLocalTraction(real strength) {
  if (local.tracEla > std::abs(strength)) {
    local.xyTraction = ((local.sXY + local.xyStress) / local.tracEla) * strength;
    local.xzTraction = ((local.sXZ + local.xzStress) / local.tracEla) * strength;

    // update stress change
    local.xyTraction -= local.sXY;
    local.xzTraction -= local.sXZ;
  } else {
    local.xyTraction = local.xyStress;
    local.xzTraction = local.xzStress;
  }
}

void ReceiverBasedOutput::computeSlipAndRate(std::array<real, 6>& rotatedTraction,
                                             std::array<real, 6>& rotatedLocalStress) {

  auto& impAndEta = ((local.layer->var(drDescr->impAndEta))[local.ltsId]);
  local.srS = -impAndEta.invEtaS * (rotatedTraction[3] - rotatedLocalStress[3]);
  local.srD = -impAndEta.invEtaS * (rotatedTraction[5] - rotatedLocalStress[5]);
}

void ReceiverBasedOutput::computeSlipAndRate(const double* tangent1,
                                             const double* tangent2,
                                             const double* strike,
                                             const double* dip) {
  local.srS = static_cast<real>(0.0);
  local.srD = static_cast<real>(0.0);

  for (size_t i = 0; i < 3; ++i) {
    real factorMinus = (local.faceAlignedValuesMinus[7] * tangent1[i] +
                        local.faceAlignedValuesMinus[8] * tangent2[i]);

    real factorPlus = (local.faceAlignedValuesPlus[7] * tangent1[i] +
                       local.faceAlignedValuesPlus[8] * tangent2[i]);

    local.srS += factorMinus * strike[i] - factorPlus * strike[i];
    local.srD += factorMinus * dip[i] - factorPlus * dip[i];
  }
}

int ReceiverBasedOutput::getClosestInternalGp(int nearestGpIndex, int nPoly) {
  int i1 = int((nearestGpIndex - 1) / (nPoly + 2)) + 1;
  int j1 = (nearestGpIndex - 1) % (nPoly + 2) + 1;
  if (i1 == 1) {
    i1 = i1 + 1;
  } else if (i1 == (nPoly + 2)) {
    i1 = i1 - 1;
  }

  if (j1 == 1) {
    j1 = j1 + 1;
  } else if (j1 == (nPoly + 2)) {
    j1 = j1 - 1;
  }
  return (i1 - 1) * (nPoly + 2) + j1;
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
      computeTriDubinerPolynomials(phiAtPoint.data(), chi, tau, numPoly);

      for (size_t d = 0; d < numDegFr2d; ++d) {
        projectedRT[d] += weights(jBndGP) * rt[local.ltsId][jBndGP] * phiAtPoint[d];
      }
    }

    auto m2inv =
        seissol::init::M2inv::view::create(const_cast<real*>(seissol::init::M2inv::Values));
    for (size_t d = 0; d < numDegFr2d; ++d) {
      projectedRT[d] *= m2inv(d, d);
    }

    auto nearestGpIndex = getClosestInternalGp(local.nearestGpIndex, numPoly);
    real chi = chiTau2dPoints(nearestGpIndex, 0);
    real tau = chiTau2dPoints(nearestGpIndex, 1);

    computeGradTriDubinerPolynomials(phiAtPoint.data(), chi, tau, numPoly);

    real dTdChi{0.0};
    real dTdTau{0.0};
    for (size_t d = 0; d < numDegFr2d; ++d) {
      dTdChi += projectedRT[d] * phiAtPoint[2 * d];
      dTdTau += projectedRT[d] * phiAtPoint[2 * d + 1];
    }
    real dTdX = jacobiT2d(0, 0) * dTdChi + jacobiT2d(0, 1) * dTdTau;
    real dTdY = jacobiT2d(1, 0) * dTdChi + jacobiT2d(1, 1) * dTdTau;

    real slowness = std::sqrt(dTdX * dTdX + dTdY * dTdY);
    ruptureVelocity = (slowness == 0.0) ? 0.0 : 1.0 / slowness;
  }

  return ruptureVelocity;
}

void ReceiverBasedOutput::tiePointers(seissol::initializers::Layer& layerData,
                                      seissol::initializers::DynamicRupture* description,
                                      seissol::Interoperability& eInteroperability) {
  constexpr auto size = init::QInterpolated::Stop[0];
  real(*accumulatedSlipMagnitude)[size] = layerData.var(description->accumulatedSlipMagnitude);
  real(*slip1)[size] = layerData.var(description->slip1);
  real(*slip2)[size] = layerData.var(description->slip2);
  real(*ruptureTime)[size] = layerData.var(description->ruptureTime);
  real(*dynStressTime)[size] = layerData.var(description->dynStressTime);
  real(*peakSR)[size] = layerData.var(description->peakSlipRate);
  real(*traction1)[size] = layerData.var(description->traction1);
  real(*traction2)[size] = layerData.var(description->traction2);

  DRFaceInformation* faceInformation = layerData.var(description->faceInformation);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
    unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
    eInteroperability.copyFrictionOutputToFortranGeneral(ltsFace,
                                                         meshFace,
                                                         accumulatedSlipMagnitude,
                                                         slip1,
                                                         slip2,
                                                         ruptureTime,
                                                         dynStressTime,
                                                         peakSR,
                                                         traction1,
                                                         traction2);
  }
}
} // namespace seissol::dr::output
