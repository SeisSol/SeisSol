#include "Initializer/tree/Layer.hpp"
#include "Initializer/preProcessorMacros.hpp"
#include "Numerical_aux/BasisFunction.h"
#include "ReceiverBasedOutput.hpp"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"
#include <cstring>

using namespace seissol::dr::misc::quantity_indices;

namespace seissol::dr::output {
void ReceiverOutput::setLtsData(seissol::initializers::LTSTree* userWpTree,
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

void ReceiverOutput::getDofs(real dofs[tensor::Q::size()], int meshId) {
  // get DOFs from 0th derivatives
  assert((wpLut->lookup(wpDescr->cellInformation, meshId).ltsSetup >> 9) % 2 == 1);

  real* derivatives = wpLut->lookup(wpDescr->derivatives, meshId);
#ifdef ACL_DEVICE
  device::DeviceInstance::getInstance().api->copyFrom(
      &dofs[0], &derivatives[0], sizeof(real) * tensor::dQ::Size[0]);
#else
  std::copy(&derivatives[0], &derivatives[tensor::dQ::Size[0]], &dofs[0]);
#endif
}

void ReceiverOutput::getNeighbourDofs(real dofs[tensor::Q::size()], int meshId, int side) {
  real* derivatives = wpLut->lookup(wpDescr->faceNeighbors, meshId)[side];
  assert(derivatives != nullptr);

#ifdef ACL_DEVICE
  device::DeviceInstance::getInstance().api->copyFrom(
      &dofs[0], &derivatives[0], sizeof(real) * tensor::dQ::Size[0]);
#else
  std::copy(&derivatives[0], &derivatives[tensor::dQ::Size[0]], &dofs[0]);
#endif
}

void ReceiverOutput::allocateMemory(
    const std::vector<std::shared_ptr<ReceiverOutputData>>& states) {
#ifdef ACL_DEVICE
  std::size_t maxCellCount = 0;
  for (const auto& state : states) {
    if (state) {
      maxCellCount = std::max(state->cellCount, maxCellCount);
    }
  }
  deviceCopyMemory =
      reinterpret_cast<real*>(device::DeviceInstance::getInstance().api->allocPinnedMem(
          sizeof(real) * tensor::Q::size() * maxCellCount));
#endif
}

ReceiverOutput::~ReceiverOutput() {
  if (deviceCopyMemory != nullptr) {
#ifdef ACL_DEVICE
    device::DeviceInstance::getInstance().api->freePinnedMem(deviceCopyMemory);
#endif
    deviceCopyMemory = nullptr;
  }
}

void ReceiverOutput::calcFaultOutput(const OutputType type,
                                     std::shared_ptr<ReceiverOutputData> outputData,
                                     const GeneralParams& generalParams,
                                     double time) {

  const size_t level = (type == OutputType::AtPickpoint) ? outputData->currentCacheLevel : 0;
  const auto faultInfos = meshReader->getFault();

#ifdef ACL_DEVICE
  if (outputData->cellCount > 0) {
    void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();
    device::DeviceInstance::getInstance().algorithms.copyScatterToUniform(outputData->deviceDataPtr,
                                                                          deviceCopyMemory,
                                                                          tensor::Q::size(),
                                                                          tensor::Q::size(),
                                                                          outputData->cellCount,
                                                                          stream);
    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
  }
#endif

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < outputData->receiverPoints.size(); ++i) {
    alignas(ALIGNMENT) real dofsPlus[tensor::Q::size()]{};
    alignas(ALIGNMENT) real dofsMinus[tensor::Q::size()]{};

    assert(outputData->receiverPoints[i].isInside == true &&
           "a receiver is not within any tetrahedron adjacent to a fault");

    const auto faceIndex = outputData->receiverPoints[i].faultFaceIndex;
    assert(faceIndex != -1 && "receiver is not initialized");
    LocalInfo local{};

    auto [layer, ltsId] = (*faceToLtsMap)[faceIndex];
    local.layer = layer;
    local.ltsId = ltsId;

    local.nearestGpIndex = outputData->receiverPoints[i].nearestGpIndex;
    local.nearestInternalGpIndex = outputData->receiverPoints[i].nearestInternalGpIndex;

    local.waveSpeedsPlus = &((local.layer->var(drDescr->waveSpeedsPlus))[local.ltsId]);
    local.waveSpeedsMinus = &((local.layer->var(drDescr->waveSpeedsMinus))[local.ltsId]);

    const auto faultInfo = faultInfos[faceIndex];

#ifdef ACL_DEVICE
    {
      real* dofsPlusData = deviceCopyMemory + tensor::Q::size() * outputData->deviceDataPlus[i];
      real* dofsMinusData = deviceCopyMemory + tensor::Q::size() * outputData->deviceDataMinus[i];

      std::memcpy(dofsPlus, dofsPlusData, sizeof(dofsPlus));
      std::memcpy(dofsMinus, dofsMinusData, sizeof(dofsMinus));
    }
#else
    getDofs(dofsPlus, faultInfo.element);
    if (faultInfo.neighborElement >= 0) {
      getDofs(dofsMinus, faultInfo.neighborElement);
    } else {
      getNeighbourDofs(dofsMinus, faultInfo.element, faultInfo.side);
    }
#endif

    const auto* initStresses = local.layer->var(drDescr->initialStressInFaultCS);
    const auto* initStress = initStresses[local.ltsId][local.nearestGpIndex];

    local.frictionCoefficient = (local.layer->var(drDescr->mu))[local.ltsId][local.nearestGpIndex];
    local.stateVariable = this->computeStateVariable(local);

    local.iniTraction1 = initStress[QuantityIndices::XY];
    local.iniTraction2 = initStress[QuantityIndices::XZ];
    local.iniNormalTraction = initStress[QuantityIndices::XX];
    local.fluidPressure = this->computeFluidPressure(local);

    const auto& normal = outputData->faultDirections[i].faceNormal;
    const auto& tangent1 = outputData->faultDirections[i].tangent1;
    const auto& tangent2 = outputData->faultDirections[i].tangent2;
    const auto& strike = outputData->faultDirections[i].strike;
    const auto& dip = outputData->faultDirections[i].dip;

    auto* phiPlusSide = outputData->basisFunctions[i].plusSide.data();
    auto* phiMinusSide = outputData->basisFunctions[i].minusSide.data();

    seissol::dynamicRupture::kernel::evaluateFaceAlignedDOFSAtPoint kernel;
    kernel.Tinv = outputData->glbToFaceAlignedData[i].data();

    kernel.Q = dofsPlus;
    kernel.basisFunctionsAtPoint = phiPlusSide;
    kernel.QAtPoint = local.faceAlignedValuesPlus;
    kernel.execute();

    kernel.Q = dofsMinus;
    kernel.basisFunctionsAtPoint = phiMinusSide;
    kernel.QAtPoint = local.faceAlignedValuesMinus;
    kernel.execute();

    this->computeLocalStresses(local);
    const real strength = this->computeLocalStrength(local);
    this->updateLocalTractions(local, strength);

    seissol::dynamicRupture::kernel::rotateInitStress alignAlongDipAndStrikeKernel;
    alignAlongDipAndStrikeKernel.stressRotationMatrix =
        outputData->stressGlbToDipStrikeAligned[i].data();
    alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
        outputData->stressFaceAlignedToGlb[i].data();

    std::array<real, 6> updatedStress{};
    updatedStress[QuantityIndices::XX] = local.transientNormalTraction;
    updatedStress[QuantityIndices::YY] = local.faceAlignedStress22;
    updatedStress[QuantityIndices::ZZ] = local.faceAlignedStress33;
    updatedStress[QuantityIndices::XY] = local.updatedTraction1;
    updatedStress[QuantityIndices::YZ] = local.faceAlignedStress23;
    updatedStress[QuantityIndices::XZ] = local.updatedTraction2;

    alignAlongDipAndStrikeKernel.initialStress = updatedStress.data();
    std::array<real, 6> rotatedUpdatedStress{};
    alignAlongDipAndStrikeKernel.rotatedStress = rotatedUpdatedStress.data();
    alignAlongDipAndStrikeKernel.execute();

    std::array<real, 6> stress{};
    stress[QuantityIndices::XX] = local.transientNormalTraction;
    stress[QuantityIndices::YY] = local.faceAlignedStress22;
    stress[QuantityIndices::ZZ] = local.faceAlignedStress33;
    stress[QuantityIndices::XY] = local.faceAlignedStress12;
    stress[QuantityIndices::YZ] = local.faceAlignedStress23;
    stress[QuantityIndices::XZ] = local.faceAlignedStress13;

    alignAlongDipAndStrikeKernel.initialStress = stress.data();
    std::array<real, 6> rotatedStress{};
    alignAlongDipAndStrikeKernel.rotatedStress = rotatedStress.data();
    alignAlongDipAndStrikeKernel.execute();

    switch (generalParams.slipRateOutputType) {
    case SlipRateOutputType::TractionsAndFailure: {
      this->computeSlipRate(local, rotatedUpdatedStress, rotatedStress);
      break;
    }
    case SlipRateOutputType::VelocityDifference: {
      this->computeSlipRate(local, tangent1, tangent2, strike, dip);
      break;
    }
    }

    adjustRotatedUpdatedStress(rotatedUpdatedStress, rotatedStress);

    auto& slipRate = std::get<VariableID::SlipRate>(outputData->vars);
    if (slipRate.isActive) {
      slipRate(DirectionID::Strike, level, i) = local.slipRateStrike;
      slipRate(DirectionID::Dip, level, i) = local.slipRateDip;
    }

    auto& transientTractions = std::get<VariableID::TransientTractions>(outputData->vars);
    if (transientTractions.isActive) {
      transientTractions(DirectionID::Strike, level, i) = rotatedUpdatedStress[QuantityIndices::XY];
      transientTractions(DirectionID::Dip, level, i) = rotatedUpdatedStress[QuantityIndices::XZ];
      transientTractions(DirectionID::Normal, level, i) =
          local.transientNormalTraction - local.fluidPressure;
    }

    auto& frictionAndState = std::get<VariableID::FrictionAndState>(outputData->vars);
    if (frictionAndState.isActive) {
      frictionAndState(ParamID::FrictionCoefficient, level, i) = local.frictionCoefficient;
      frictionAndState(ParamID::State, level, i) = local.stateVariable;
    }

    auto& ruptureTime = std::get<VariableID::RuptureTime>(outputData->vars);
    if (ruptureTime.isActive) {
      auto* rt = local.layer->var(drDescr->ruptureTime);
      ruptureTime(level, i) = rt[local.ltsId][local.nearestGpIndex];
    }

    auto& normalVelocity = std::get<VariableID::NormalVelocity>(outputData->vars);
    if (normalVelocity.isActive) {
      normalVelocity(level, i) = local.faultNormalVelocity;
    }

    auto& accumulatedSlip = std::get<VariableID::AccumulatedSlip>(outputData->vars);
    if (accumulatedSlip.isActive) {
      auto* slip = local.layer->var(drDescr->accumulatedSlipMagnitude);
      accumulatedSlip(level, i) = slip[local.ltsId][local.nearestGpIndex];
    }

    auto& totalTractions = std::get<VariableID::TotalTractions>(outputData->vars);
    if (totalTractions.isActive) {
      std::array<real, tensor::rotatedStress::size()> rotatedInitStress{};
      alignAlongDipAndStrikeKernel.initialStress = initStress;
      alignAlongDipAndStrikeKernel.rotatedStress = rotatedInitStress.data();
      alignAlongDipAndStrikeKernel.execute();

      totalTractions(DirectionID::Strike, level, i) =
          rotatedUpdatedStress[QuantityIndices::XY] + rotatedInitStress[QuantityIndices::XY];
      totalTractions(DirectionID::Dip, level, i) =
          rotatedUpdatedStress[QuantityIndices::XZ] + rotatedInitStress[QuantityIndices::XZ];
      totalTractions(DirectionID::Normal, level, i) = local.transientNormalTraction -
                                                      local.fluidPressure +
                                                      rotatedInitStress[QuantityIndices::XX];
    }

    auto& ruptureVelocity = std::get<VariableID::RuptureVelocity>(outputData->vars);
    if (ruptureVelocity.isActive) {
      auto& jacobiT2d = outputData->jacobianT2d[i];
      ruptureVelocity(level, i) = this->computeRuptureVelocity(jacobiT2d, local);
    }

    auto& peakSlipsRate = std::get<VariableID::PeakSlipRate>(outputData->vars);
    if (peakSlipsRate.isActive) {
      auto* peakSR = local.layer->var(drDescr->peakSlipRate);
      peakSlipsRate(level, i) = peakSR[local.ltsId][local.nearestGpIndex];
    }

    auto& dynamicStressTime = std::get<VariableID::DynamicStressTime>(outputData->vars);
    if (dynamicStressTime.isActive) {
      auto* dynStressTime = (local.layer->var(drDescr->dynStressTime));
      dynamicStressTime(level, i) = dynStressTime[local.ltsId][local.nearestGpIndex];
    }

    auto& slipVectors = std::get<VariableID::Slip>(outputData->vars);
    if (slipVectors.isActive) {
      VrtxCoords crossProduct = {0.0, 0.0, 0.0};
      MeshTools::cross(strike.data(), tangent1.data(), crossProduct);

      const double cos1 = MeshTools::dot(strike.data(), tangent1.data());
      const double scalarProd = MeshTools::dot(crossProduct, normal.data());

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
    this->outputSpecifics(outputData, local, level, i);
  }

  if (type == OutputType::AtPickpoint) {
    outputData->cachedTime[outputData->currentCacheLevel] = time;
    outputData->currentCacheLevel += 1;
  }
}

void ReceiverOutput::computeLocalStresses(LocalInfo& local) {
  const auto& impAndEta = ((local.layer->var(drDescr->impAndEta))[local.ltsId]);
  const real normalDivisor = 1.0 / (impAndEta.zpNeig + impAndEta.zp);
  const real shearDivisor = 1.0 / (impAndEta.zsNeig + impAndEta.zs);

  auto diff = [&local](int i) {
    return local.faceAlignedValuesMinus[i] - local.faceAlignedValuesPlus[i];
  };

  local.faceAlignedStress12 =
      local.faceAlignedValuesPlus[QuantityIndices::XY] +
      ((diff(QuantityIndices::XY) + impAndEta.zsNeig * diff(QuantityIndices::V)) * impAndEta.zs) *
          shearDivisor;

  local.faceAlignedStress13 =
      local.faceAlignedValuesPlus[QuantityIndices::XZ] +
      ((diff(QuantityIndices::XZ) + impAndEta.zsNeig * diff(QuantityIndices::W)) * impAndEta.zs) *
          shearDivisor;

  local.transientNormalTraction =
      local.faceAlignedValuesPlus[QuantityIndices::XX] +
      ((diff(QuantityIndices::XX) + impAndEta.zpNeig * diff(QuantityIndices::U)) * impAndEta.zp) *
          normalDivisor;

  local.faultNormalVelocity =
      local.faceAlignedValuesPlus[QuantityIndices::U] +
      (local.transientNormalTraction - local.faceAlignedValuesPlus[QuantityIndices::XX]) *
          impAndEta.invZp;

  real missingSigmaValues =
      (local.transientNormalTraction - local.faceAlignedValuesPlus[QuantityIndices::XX]);
  missingSigmaValues *= (1.0 - 2.0 * std::pow(local.waveSpeedsPlus->sWaveVelocity /
                                                  local.waveSpeedsPlus->pWaveVelocity,
                                              2));

  local.faceAlignedStress22 = local.faceAlignedValuesPlus[QuantityIndices::YY] + missingSigmaValues;
  local.faceAlignedStress33 = local.faceAlignedValuesPlus[QuantityIndices::ZZ] + missingSigmaValues;
  local.faceAlignedStress23 = local.faceAlignedValuesPlus[QuantityIndices::YZ];
}

void ReceiverOutput::updateLocalTractions(LocalInfo& local, real strength) {
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

void ReceiverOutput::computeSlipRate(LocalInfo& local,
                                     const std::array<real, 6>& rotatedUpdatedStress,
                                     const std::array<real, 6>& rotatedStress) {

  const auto& impAndEta = ((local.layer->var(drDescr->impAndEta))[local.ltsId]);
  local.slipRateStrike = -impAndEta.invEtaS * (rotatedUpdatedStress[QuantityIndices::XY] -
                                               rotatedStress[QuantityIndices::XY]);
  local.slipRateDip = -impAndEta.invEtaS * (rotatedUpdatedStress[QuantityIndices::XZ] -
                                            rotatedStress[QuantityIndices::XZ]);
}

void ReceiverOutput::computeSlipRate(LocalInfo& local,
                                     const std::array<double, 3>& tangent1,
                                     const std::array<double, 3>& tangent2,
                                     const std::array<double, 3>& strike,
                                     const std::array<double, 3>& dip) {
  local.slipRateStrike = static_cast<real>(0.0);
  local.slipRateDip = static_cast<real>(0.0);

  for (size_t i = 0; i < 3; ++i) {
    const real factorMinus = (local.faceAlignedValuesMinus[QuantityIndices::V] * tangent1[i] +
                              local.faceAlignedValuesMinus[QuantityIndices::W] * tangent2[i]);

    const real factorPlus = (local.faceAlignedValuesPlus[QuantityIndices::V] * tangent1[i] +
                             local.faceAlignedValuesPlus[QuantityIndices::W] * tangent2[i]);

    local.slipRateStrike += (factorMinus - factorPlus) * strike[i];
    local.slipRateDip += (factorMinus - factorPlus) * dip[i];
  }
}

real ReceiverOutput::computeRuptureVelocity(Eigen::Matrix<real, 2, 2>& jacobiT2d,
                                            const LocalInfo& local) {
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
      basisFunction::tri_dubiner::evaluatePolynomials(phiAtPoint.data(), chi, tau, numPoly);

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

    basisFunction::tri_dubiner::evaluateGradPolynomials(phiAtPoint.data(), chi, tau, numPoly);

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
