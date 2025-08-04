// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ReceiverBasedOutput.h"
#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Numerical/BasisFunction.h"
#include <Alignment.h>
#include <GeneratedCode/init.h>
#include <Kernels/Common.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/MultipleSimulations.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <vector>

using namespace seissol::dr::misc::quantity_indices;

namespace seissol::dr::output {
void ReceiverOutput::setLtsData(LTS::Storage& userWpStorage,
                                LTS::Backmap& userWpBackmap,
                                DynamicRupture::Storage& userDrStorage) {
  wpStorage = &userWpStorage;
  wpBackmap = &userWpBackmap;
  drStorage = &userDrStorage;
}

template<typename Derived>
template<typename Cfg>
void ReceiverOutputImpl<Derived>::getDofs(Real<Cfg> dofs[tensor::Q<Cfg>::size()], int meshId) {
  const auto position = wpBackmap->get(meshId);
  auto& layer = wpStorage->layer(position.color);
  // get DOFs from 0th derivatives
  assert(layer.var<LTS::CellInformation>()[position.cell].ltsSetup.hasDerivatives());

  Real<Cfg>* derivatives = layer.var<LTS::Derivatives>(Cfg())[position.cell];
  std::copy(&derivatives[0], &derivatives[tensor::dQ<Cfg>::Size[0]], &dofs[0]);
}

template<typename Derived>
template<typename Cfg>
void ReceiverOutputImpl<Derived>::getNeighborDofs(Real<Cfg> dofs[tensor::Q<Cfg>::size()], int meshId, int side) {
  const auto position = wpBackmap->get(meshId);
  auto& layer = wpStorage->layer(position.color);
  auto* derivatives = reinterpret_cast<Real<Cfg>*>(layer.var<LTS::FaceNeighbors>()[position.cell][side]);
  assert(derivatives != nullptr);

  std::copy(&derivatives[0], &derivatives[tensor::dQ<Cfg>::Size[0]], &dofs[0]);
}

template<typename Derived>
void ReceiverOutputImpl<Derived>::calcFaultOutput(
    seissol::initializer::parameters::OutputType outputType,
    seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
    const std::shared_ptr<ReceiverOutputData>& outputData,
    parallel::runtime::StreamRuntime& runtime,
    double time) {

  const size_t level = (outputType == seissol::initializer::parameters::OutputType::AtPickpoint)
                           ? outputData->currentCacheLevel
                           : 0;
  const auto& faultInfos = meshReader->getFault();

  auto& callRuntime =
      outputData->extraRuntime.has_value() ? outputData->extraRuntime.value() : runtime;

  if constexpr (isDeviceOn()) {
    if (outputData->extraRuntime.has_value()) {
      runtime.eventSync(outputData->extraRuntime->eventRecord());
    }
    outputData->deviceDataCollector->gatherToHost(runtime.stream());
    for (auto& [_, dataCollector] : outputData->deviceVariables) {
      dataCollector->gatherToHost(runtime.stream());
    }
    if (outputData->extraRuntime.has_value()) {
      outputData->extraRuntime->eventSync(runtime.eventRecord());
    }
  }

  const auto points = outputData->receiverPoints.size();
  const auto handler = [this, outputData, &faultInfos, outputType, slipRateOutputType, level](
                           std::size_t i) {
    // TODO: query the dofs, only once per simulation; once per face

    assert(outputData->receiverPoints[i].isInside == true &&
           "a receiver is not within any tetrahedron adjacent to a fault");

    const auto faceIndex = outputData->receiverPoints[i].faultFaceIndex;
    assert(faceIndex != -1 && "receiver is not initialized");

    auto [layer, ltsId] = (*faceToLtsMap)[faceIndex];
    layer->wrap([&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>;

    LocalInfo local{};

    alignas(Alignment) real dofsPlus[tensor::Q<Cfg>::size()]{};
    alignas(Alignment) real dofsMinus[tensor::Q<Cfg>::size()]{};

    local.layer = layer;
    local.ltsId = ltsId;
    local.index = i;
    local.fusedIndex = outputData->receiverPoints[i].simIndex;
    local.state = outputData.get();

    local.nearestGpIndex = outputData->receiverPoints[i].nearestGpIndex;
    local.gpIndex = outputData->receiverPoints[i].gpIndex;
    local.nearestInternalGpIndex = outputData->receiverPoints[i].nearestInternalGpIndex;
    local.internalGpIndexFused = outputData->receiverPoints[i].internalGpIndexFused;

    local.waveSpeedsPlus = &((local.layer->template var<DynamicRupture::WaveSpeedsPlus>())[local.ltsId]);
    local.waveSpeedsMinus = &((local.layer->template var<DynamicRupture::WaveSpeedsMinus>())[local.ltsId]);

    const auto& faultInfo = faultInfos[faceIndex];

    if constexpr (isDeviceOn()) {
      real* dofsPlusData = outputData->deviceDataCollector->get(outputData->deviceDataPlus[i]);
      real* dofsMinusData = outputData->deviceDataCollector->get(outputData->deviceDataMinus[i]);

      std::memcpy(dofsPlus, dofsPlusData, sizeof(dofsPlus));
      std::memcpy(dofsMinus, dofsMinusData, sizeof(dofsMinus));
    } else {
      getDofs(dofsPlus, faultInfo.element);
      if (faultInfo.neighborElement >= 0) {
        getDofs(dofsMinus, faultInfo.neighborElement);
      } else {
        getNeighborDofs(dofsMinus, faultInfo.element, faultInfo.side);
      }
    }

    const auto* initStresses = getCellData<DynamicRupture::InitialStressInFaultCS>(cfg, local);

    local.frictionCoefficient = getCellData<DynamicRupture::Mu>(cfg, local)[local.gpIndex];
    local.stateVariable = Derived::template computeStateVariable<Cfg>(local);

    local.iniTraction1 = initStresses[QuantityIndices::XY][local.gpIndex];
    local.iniTraction2 = initStresses[QuantityIndices::XZ][local.gpIndex];
    local.iniNormalTraction = initStresses[QuantityIndices::XX][local.gpIndex];
    local.fluidPressure = Derived::template computeFluidPressure<Cfg>(local);

    const auto& normal = outputData->faultDirections[i].faceNormal;
    const auto& tangent1 = outputData->faultDirections[i].tangent1;
    const auto& tangent2 = outputData->faultDirections[i].tangent2;
    const auto& strike = outputData->faultDirections[i].strike;
    const auto& dip = outputData->faultDirections[i].dip;

    const auto& transformData = std::get<TransformData<Cfg>>(outputData->transformData);

    auto* phiPlusSide = transformData.basisFunctions[i].plusSide.data();
    auto* phiMinusSide = transformData.basisFunctions[i].minusSide.data();

    seissol::dynamicRupture::kernel::evaluateFaceAlignedDOFSAtPoint<Cfg> kernel;
    kernel.Tinv = transformData.glbToFaceAlignedData[i].data();

    real faceAlignedValuesPlus[tensor::QAtPoint<Cfg>::size()]{};
    real faceAlignedValuesMinus[tensor::QAtPoint<Cfg>::size()]{};

    // TODO: do these operations only once per simulation
    kernel.Q = dofsPlus;
    kernel.basisFunctionsAtPoint = phiPlusSide;
    kernel.QAtPoint = faceAlignedValuesPlus;
    kernel.execute();

    kernel.Q = dofsMinus;
    kernel.basisFunctionsAtPoint = phiMinusSide;
    kernel.QAtPoint = faceAlignedValuesMinus;
    kernel.execute();

    for (size_t j = 0; j < tensor::QAtPoint<Cfg>::Shape[seissol::multisim::BasisFunctionDimension];
         ++j) {
      local.faceAlignedValuesPlus[j] =
          faceAlignedValuesPlus[j * seissol::multisim::NumSimulations + local.fusedIndex];
      local.faceAlignedValuesMinus[j] =
          faceAlignedValuesMinus[j * seissol::multisim::NumSimulations + local.fusedIndex];
    }

    this->computeLocalStresses<Cfg>(local);
    const real strength = Derived::template computeLocalStrength<Cfg>(local);
    seissol::dr::output::ReceiverOutputImpl<Derived>::updateLocalTractions<Cfg>(local, strength);

    seissol::dynamicRupture::kernel::rotateInitStress<Cfg> alignAlongDipAndStrikeKernel;
    alignAlongDipAndStrikeKernel.stressRotationMatrix =
        transformData.stressGlbToDipStrikeAligned[i].data();
    alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
        transformData.stressFaceAlignedToGlb[i].data();

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

    switch (slipRateOutputType) {
    case seissol::initializer::parameters::SlipRateOutputType::TractionsAndFailure: {
      this->computeSlipRate<Cfg>(local, rotatedUpdatedStress, rotatedStress);
      break;
    }
    case seissol::initializer::parameters::SlipRateOutputType::VelocityDifference: {
      this->computeSlipRate<Cfg>(local, tangent1, tangent2, strike, dip);
      break;
    }
    }

    Derived::template adjustRotatedUpdatedStress<Cfg>(rotatedUpdatedStress, rotatedStress);

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
      auto* rt = getCellData<DynamicRupture::RuptureTime>(cfg, local);
      ruptureTime(level, i) = rt[local.gpIndex];
    }

    auto& normalVelocity = std::get<VariableID::NormalVelocity>(outputData->vars);
    if (normalVelocity.isActive) {
      normalVelocity(level, i) = local.faultNormalVelocity;
    }

    auto& accumulatedSlip = std::get<VariableID::AccumulatedSlip>(outputData->vars);
    if (accumulatedSlip.isActive) {
      auto* slip = getCellData<DynamicRupture::AccumulatedSlipMagnitude>(cfg, local);
      accumulatedSlip(level, i) = slip[local.gpIndex];
    }

    auto& totalTractions = std::get<VariableID::TotalTractions>(outputData->vars);
    if (totalTractions.isActive) {
      std::array<real, tensor::initialStress<Cfg>::size()> unrotatedInitStress{};
      std::array<real, tensor::rotatedStress<Cfg>::size()> rotatedInitStress{};
      for (std::size_t stressVar = 0; stressVar < unrotatedInitStress.size(); ++stressVar) {
        unrotatedInitStress[stressVar] = initStresses[stressVar][local.gpIndex];
      }
      alignAlongDipAndStrikeKernel.initialStress = unrotatedInitStress.data();
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
      auto& jacobiT2d = transformData.jacobianT2d[i];
      ruptureVelocity(level, i) = this->computeRuptureVelocity<Cfg>(jacobiT2d, local);
    }

    auto& peakSlipsRate = std::get<VariableID::PeakSlipRate>(outputData->vars);
    if (peakSlipsRate.isActive) {
      auto* peakSR = getCellData<DynamicRupture::PeakSlipRate>(cfg, local);
      peakSlipsRate(level, i) = peakSR[local.gpIndex];
    }

    auto& dynamicStressTime = std::get<VariableID::DynamicStressTime>(outputData->vars);
    if (dynamicStressTime.isActive) {
      auto* dynStressTime = getCellData<DynamicRupture::DynStressTime>(cfg, local);
      dynamicStressTime(level, i) = dynStressTime[local.gpIndex];
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

      auto* slip1 = getCellData<DynamicRupture::Slip1>(cfg, local);
      auto* slip2 = getCellData<DynamicRupture::Slip2>(cfg, local);

      slipVectors(DirectionID::Strike, level, i) =
          cos1 * slip1[local.gpIndex] - sin1 * slip2[local.gpIndex];

      slipVectors(DirectionID::Dip, level, i) =
          sin1 * slip1[local.gpIndex] + cos1 * slip2[local.gpIndex];
    }
    Derived::template outputSpecifics<Cfg>(outputData, local, level, i);
  });
  };

  callRuntime.enqueueLoop(points, handler);

  if (outputType == seissol::initializer::parameters::OutputType::AtPickpoint) {
    outputData->cachedTime[outputData->currentCacheLevel] = time;
    outputData->currentCacheLevel += 1;
  }
}

template<typename Derived>
template<typename Cfg>
void ReceiverOutputImpl<Derived>::computeLocalStresses(LocalInfo<Cfg>& local) {
  using real = Real<Cfg>;

  const auto& impAndEta = ((local.layer->template var<DynamicRupture::ImpAndEta>(Cfg()))[local.ltsId]);
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

template<typename Derived>
template<typename Cfg>
void ReceiverOutputImpl<Derived>::updateLocalTractions(LocalInfo<Cfg>& local, Real<Cfg> strength) {
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

template<typename Derived>
template<typename Cfg>
void ReceiverOutputImpl<Derived>::computeSlipRate(LocalInfo<Cfg>& local,
                                     const std::array<Real<Cfg>, 6>& rotatedUpdatedStress,
                                     const std::array<Real<Cfg>, 6>& rotatedStress) {

  const auto& impAndEta = ((local.layer->template var<DynamicRupture::ImpAndEta>(Cfg()))[local.ltsId]);
  local.slipRateStrike = -impAndEta.invEtaS * (rotatedUpdatedStress[QuantityIndices::XY] -
                                               rotatedStress[QuantityIndices::XY]);
  local.slipRateDip = -impAndEta.invEtaS * (rotatedUpdatedStress[QuantityIndices::XZ] -
                                            rotatedStress[QuantityIndices::XZ]);
}

template<typename Derived>
template<typename Cfg>
void ReceiverOutputImpl<Derived>::computeSlipRate(LocalInfo<Cfg>& local,
                                     const std::array<double, 3>& tangent1,
                                     const std::array<double, 3>& tangent2,
                                     const std::array<double, 3>& strike,
                                     const std::array<double, 3>& dip) {
  using real = Real<Cfg>;
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

template<typename Derived>
template<typename Cfg>
Real<Cfg> ReceiverOutputImpl<Derived>::computeRuptureVelocity(Eigen::Matrix<Real<Cfg>, 2, 2>& jacobiT2d,
                                            const LocalInfo<Cfg>& local) {
  using real = Real<Cfg>;
  auto* ruptureTime = getCellData<DynamicRupture::RuptureTime>(Cfg(), local);
  real ruptureVelocity = 0.0;

  bool needsUpdate{true};
  for (size_t point = 0; point < misc::NumBoundaryGaussPoints<Cfg>; ++point) {
    if (ruptureTime[point] == 0.0) {
      needsUpdate = false;
    }
  }

  if (needsUpdate) {
    constexpr int NumPoly = Cfg::ConvergenceOrder - 1;
    constexpr int NumDegFr2d = (NumPoly + 1) * (NumPoly + 2) / 2;
    std::array<double, NumDegFr2d> projectedRT{};
    projectedRT.fill(0.0);

    std::array<double, static_cast<std::size_t>(2 * NumDegFr2d)> phiAtPoint{};
    phiAtPoint.fill(0.0);

    auto chiTau2dPoints =
        init::quadpoints<Cfg>::view::create(const_cast<real*>(init::quadpoints<Cfg>::Values));
    auto weights =
        init::quadweights<Cfg>::view::create(const_cast<real*>(init::quadweights<Cfg>::Values));

    auto* rt = getCellData<DynamicRupture::RuptureTime>(Cfg(), local);
    for (size_t jBndGP = 0; jBndGP < misc::NumBoundaryGaussPoints<Cfg>; ++jBndGP) {
      const real chi = seissol::multisim::multisimTranspose(chiTau2dPoints, jBndGP, 0);
      const real tau = seissol::multisim::multisimTranspose(chiTau2dPoints, jBndGP, 1);

      basisFunction::tri_dubiner::evaluatePolynomials(phiAtPoint.data(), chi, tau, NumPoly);

      for (size_t d = 0; d < NumDegFr2d; ++d) {
        projectedRT[d] +=
            seissol::multisim::multisimWrap(weights, 0, jBndGP) * rt[jBndGP] * phiAtPoint[d];
      }
    }
    auto m2inv = seissol::init::M2inv<Cfg>::view::create(
        const_cast<real*>(seissol::init::M2inv<Cfg>::Values));
    for (size_t d = 0; d < NumDegFr2d; ++d) {
      projectedRT[d] *= m2inv(d, d);
    }

    const real chi =
        seissol::multisim::multisimTranspose(chiTau2dPoints, local.nearestInternalGpIndex, 0);
    const real tau =
        seissol::multisim::multisimTranspose(chiTau2dPoints, local.nearestInternalGpIndex, 1);
    basisFunction::tri_dubiner::evaluateGradPolynomials(phiAtPoint.data(), chi, tau, NumPoly);

    real dTdChi{0.0};
    real dTdTau{0.0};
    for (size_t d = 0; d < NumDegFr2d; ++d) {
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

template<typename Derived>
std::vector<std::size_t> ReceiverOutputImpl<Derived>::getOutputVariables() const {
  return {drStorage->info<DynamicRupture::InitialStressInFaultCS>().index,
          drStorage->info<DynamicRupture::Mu>().index,
          drStorage->info<DynamicRupture::RuptureTime>().index,
          drStorage->info<DynamicRupture::AccumulatedSlipMagnitude>().index,
          drStorage->info<DynamicRupture::PeakSlipRate>().index,
          drStorage->info<DynamicRupture::DynStressTime>().index,
          drStorage->info<DynamicRupture::Slip1>().index,
          drStorage->info<DynamicRupture::Slip2>().index};
}

} // namespace seissol::dr::output
