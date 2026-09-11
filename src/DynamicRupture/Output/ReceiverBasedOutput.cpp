// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ReceiverBasedOutput.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Typedefs.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Quadrature.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
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
  wpStorage_ = &userWpStorage;
  wpBackmap_ = &userWpBackmap;
  drStorage_ = &userDrStorage;
}

void ReceiverOutput::getDofs(const real*(&derivatives), std::size_t meshId) {
  const auto position = wpBackmap_->get(meshId);
  auto& layer = wpStorage_->layer(position.color);
  // get DOFs from 0th derivatives
  assert(layer.var<LTS::CellInformation>()[position.cell].ltsSetup.hasDerivatives());

  derivatives = layer.var<LTS::Derivatives>()[position.cell];
}

void ReceiverOutput::getNeighborDofs(const real*(&derivatives),
                                     std::size_t meshId,
                                     std::size_t side) {
  const auto position = wpBackmap_->get(meshId);
  auto& layer = wpStorage_->layer(position.color);

  derivatives = layer.var<LTS::FaceNeighbors>()[position.cell][side];
  assert(derivatives != nullptr);
}

void ReceiverOutput::calcFaultOutput(
    seissol::initializer::parameters::OutputType outputType,
    seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
    const std::shared_ptr<ReceiverOutputData>& outputData,
    parallel::runtime::StreamRuntime& runtime,
    double time,
    double dt,
    double indt) {

  const size_t level = (outputType == seissol::initializer::parameters::OutputType::AtPickpoint)
                           ? outputData->currentCacheLevel
                           : 0;
  const auto& faultInfos = meshReader_->getFault();

  // the friction solve advances in the sub intervals of the time quadrature; the stored friction
  // state belongs to the last of them
  const auto frictionTime = seissol::dr::friction_law::FrictionSolver::computeDeltaT(
      seissol::quadrature::ShiftedGaussLegendre(ConvergenceOrder, 0, dt).first);

  const auto timeCoeffs = kernels::timeBasis().point(indt, dt);
  auto integrateCoeffs = kernels::timeBasis().integrate(0, indt, dt);
  for (auto& coeff : integrateCoeffs) {
    coeff = -coeff;
  }

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
  const auto handler = [this,
                        outputData,
                        &faultInfos,
                        outputType,
                        slipRateOutputType,
                        level,
                        timeCoeffs,
                        integrateCoeffs,
                        time,
                        frictionTime](std::size_t i) {
    // TODO: query the dofs, only once per simulation; once per face
    alignas(Alignment) real dofsPlus[tensor::Q::size()]{};
    alignas(Alignment) real dofsMinus[tensor::Q::size()]{};

    assert(outputData->receiverPoints[i].isInside == true &&
           "a receiver is not within any tetrahedron adjacent to a fault");

    const auto faceIndex = outputData->receiverPoints[i].faultFaceIndex;
    assert(faceIndex != -1 && "receiver is not initialized");
    LocalInfo local{};

    const auto position = faceToLtsMap_->get(faceIndex);
    local.layer = &drStorage_->layer(position.color);
    local.ltsId = position.cell;
    local.index = i;
    local.fusedIndex = outputData->receiverPoints[i].simIndex;
    local.state = outputData.get();

    local.time = time;
    local.deltaT = frictionTime.deltaT.back();
    local.printWarning = &this->printRSFWarning_;

    local.nearestGpIndex = outputData->receiverPoints[i].nearestGpIndex;
    local.gpIndex = outputData->receiverPoints[i].gpIndex;
    local.nearestInternalGpIndex = outputData->receiverPoints[i].nearestInternalGpIndex;
    local.internalGpIndexFused = outputData->receiverPoints[i].internalGpIndexFused;

    local.waveSpeedsPlus = &((local.layer->var<DynamicRupture::WaveSpeedsPlus>())[local.ltsId]);
    local.waveSpeedsMinus = &((local.layer->var<DynamicRupture::WaveSpeedsMinus>())[local.ltsId]);

    this->handleNonConvergence(local);

    const auto& faultInfo = faultInfos[faceIndex];

    if (outputType == initializer::parameters::OutputType::Elementwise) {
      std::memcpy(dofsPlus,
                  local.layer->var<DynamicRupture::TimeDofsPlus>()[local.ltsId],
                  sizeof(dofsPlus));
      std::memcpy(dofsMinus,
                  local.layer->var<DynamicRupture::TimeDofsMinus>()[local.ltsId],
                  sizeof(dofsMinus));
    } else {
      // only interpolate for the on-fault receivers
      const real* stePlus = nullptr;
      const real* steMinus = nullptr;

      if constexpr (isDeviceOn()) {
        stePlus = outputData->deviceDataCollector->get(outputData->deviceDataPlus[i]);
        steMinus = outputData->deviceDataCollector->get(outputData->deviceDataMinus[i]);
      } else {
        getDofs(stePlus, faultInfo.element);
        if (faultInfo.neighborElement >= 0) {
          getDofs(steMinus, faultInfo.neighborElement);
        } else {
          getNeighborDofs(steMinus, faultInfo.element, faultInfo.side);
        }
      }

      timeKernel_.evaluate(timeCoeffs.data(), stePlus, dofsPlus);
      timeKernel_.evaluate(timeCoeffs.data(), steMinus, dofsMinus);
    }

    const auto* initStresses = getCellData<DynamicRupture::InitialStressInFaultCS>(local);

    local.frictionCoefficient = getCellData<DynamicRupture::Mu>(local)[local.gpIndex];
    local.stateVariable = this->computeStateVariable(local);

    local.iniTraction1 = initStresses[QuantityIndices::XY][local.gpIndex];
    local.iniTraction2 = initStresses[QuantityIndices::XZ][local.gpIndex];
    local.iniNormalTraction = initStresses[QuantityIndices::XX][local.gpIndex];
    local.fluidPressure = this->computeFluidPressure(local);

    const auto& normal = outputData->faultDirections[i].faceNormal;
    const auto& tangent1 = outputData->faultDirections[i].tangent1;
    const auto& tangent2 = outputData->faultDirections[i].tangent2;
    const auto& strike = outputData->faultDirections[i].strike;
    const auto& dip = outputData->faultDirections[i].dip;

    const auto* phiPlusSide = outputData->basisFunctions[i].plusSide.data();
    const auto* phiMinusSide = outputData->basisFunctions[i].minusSide.data();

    seissol::dynamicRupture::kernel::evaluateFaceAlignedDOFSAtPoint kernel;
    kernel.Tinv = outputData->glbToFaceAlignedData[i].data();

    real faceAlignedValuesPlus[tensor::QAtPoint::size()]{};
    real faceAlignedValuesMinus[tensor::QAtPoint::size()]{};

    // TODO: do these operations only once per simulation
    kernel.Q = dofsPlus;
    kernel.basisFunctionsAtPoint = phiPlusSide;
    kernel.QAtPoint = faceAlignedValuesPlus;
    kernel.execute();

    kernel.Q = dofsMinus;
    kernel.basisFunctionsAtPoint = phiMinusSide;
    kernel.QAtPoint = faceAlignedValuesMinus;
    kernel.execute();

    for (size_t j = 0; j < tensor::QAtPoint::Shape[seissol::multisim::BasisFunctionDimension];
         ++j) {
      local.faceAlignedValuesPlus[j] =
          faceAlignedValuesPlus[j * seissol::multisim::NumSimulations + local.fusedIndex];
      local.faceAlignedValuesMinus[j] =
          faceAlignedValuesMinus[j * seissol::multisim::NumSimulations + local.fusedIndex];
    }

    this->computeLocalStresses(local);
    const real strength = this->computeLocalStrength(local);
    const real strengthSlope = this->computeLocalStrengthSlope(local);
    seissol::dr::output::ReceiverOutput::updateLocalTractions(local, strength, strengthSlope);

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

    switch (slipRateOutputType) {
    case seissol::initializer::parameters::SlipRateOutputType::TractionsAndFailure: {
      this->computeSlipRate(
          local, rotatedUpdatedStress, rotatedStress, tangent1, tangent2, strike, dip);
      break;
    }
    case seissol::initializer::parameters::SlipRateOutputType::VelocityDifference: {
      seissol::dr::output::ReceiverOutput::computeSlipRate(local, tangent1, tangent2, strike, dip);
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
      const auto* rt = getCellData<DynamicRupture::RuptureTime>(local);
      ruptureTime(level, i) = rt[local.gpIndex];
    }

    auto& normalVelocity = std::get<VariableID::NormalVelocity>(outputData->vars);
    if (normalVelocity.isActive) {
      normalVelocity(level, i) = local.faultNormalVelocity;
    }

    auto& accumulatedSlip = std::get<VariableID::AccumulatedSlip>(outputData->vars);
    if (accumulatedSlip.isActive) {
      const auto* slip = getCellData<DynamicRupture::AccumulatedSlipMagnitude>(local);
      accumulatedSlip(level, i) = slip[local.gpIndex];
    }

    auto& totalTractions = std::get<VariableID::TotalTractions>(outputData->vars);
    if (totalTractions.isActive) {
      std::array<real, tensor::initialStress::size()> unrotatedInitStress{};
      std::array<real, tensor::rotatedStress::size()> rotatedInitStress{};
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
      const auto& jacobiT2d = outputData->jacobianT2d[i];
      ruptureVelocity(level, i) = this->computeRuptureVelocity(jacobiT2d, local);
    }

    auto& peakSlipsRate = std::get<VariableID::PeakSlipRate>(outputData->vars);
    if (peakSlipsRate.isActive) {
      const auto* peakSR = getCellData<DynamicRupture::PeakSlipRate>(local);
      peakSlipsRate(level, i) = peakSR[local.gpIndex];
    }

    auto& dynamicStressTime = std::get<VariableID::DynamicStressTime>(outputData->vars);
    if (dynamicStressTime.isActive) {
      const auto* dynStressTime = getCellData<DynamicRupture::DynStressTime>(local);
      dynamicStressTime(level, i) = dynStressTime[local.gpIndex];
    }

    auto& slipVectors = std::get<VariableID::Slip>(outputData->vars);
    if (slipVectors.isActive) {
      VrtxCoords crossProduct = {0.0, 0.0, 0.0};
      MeshTools::cross(strike.data(), tangent1.data(), crossProduct);

      const double cos1t = MeshTools::dot(strike.data(), tangent1.data());
      const double scalarProd = MeshTools::dot(crossProduct, normal.data());

      // Note: cos1t**2 can be greater than 1.0 because of rounding errors -> min
      double sin1t = std::sqrt(1.0 - std::min(1.0, cos1t * cos1t));
      sin1t = (scalarProd > 0) ? sin1t : -sin1t;

      const auto* slip1 = getCellData<DynamicRupture::Slip1>(local);
      const auto* slip2 = getCellData<DynamicRupture::Slip2>(local);

      slipVectors(DirectionID::Strike, level, i) =
          cos1t * slip1[local.gpIndex] - sin1t * slip2[local.gpIndex];

      slipVectors(DirectionID::Dip, level, i) =
          sin1t * slip1[local.gpIndex] + cos1t * slip2[local.gpIndex];
    }
    this->outputSpecifics(outputData, local, level, i);
  };

  callRuntime.enqueueLoop(points, handler);

  if (outputType == seissol::initializer::parameters::OutputType::AtPickpoint) {
    outputData->cachedTime[outputData->currentCacheLevel] = time;
    outputData->currentCacheLevel += 1;
  }
}

void ReceiverOutput::computeLocalStresses(LocalInfo& local) {
  auto diff = [&local](int i) {
    return local.faceAlignedValuesMinus[i] - local.faceAlignedValuesPlus[i];
  };

  if constexpr (model::MaterialT::Type != model::MaterialType::Elastic &&
                model::MaterialT::Type != model::MaterialType::Viscoelastic) {
    // Anisotropy couples the fault-normal and the two tangential directions, poroelasticity adds
    // the fluid pressure as a fourth interface variable. In both cases the Godunov state has to be
    // reconstructed with the full matrix -- exactly as
    // common::precomputeStressFromQInterpolated does it for the solver:
    //
    //   T*   = eta * (Y+ T+ + Y- T- + (v- - v+))
    //   v*_+ = v+ + Y+ (T* - T+)
    const auto& impedanceMatrices =
        ((local.layer->var<DynamicRupture::ImpedanceMatrices>())[local.ltsId]);

    constexpr std::size_t Count =
        model::MaterialT::Type == model::MaterialType::Poroelastic ? 4 : 3;
    constexpr auto StressIndices = []() {
      if constexpr (Count == 4) {
        return std::array<int, 4>{
            QuantityIndices::XX, QuantityIndices::XY, QuantityIndices::XZ, QuantityIndices::FP};
      } else {
        return std::array<int, 3>{QuantityIndices::XX, QuantityIndices::XY, QuantityIndices::XZ};
      }
    }();
    constexpr auto VelocityIndices = []() {
      if constexpr (Count == 4) {
        return std::array<int, 4>{
            QuantityIndices::U, QuantityIndices::V, QuantityIndices::W, QuantityIndices::FU};
      } else {
        return std::array<int, 3>{QuantityIndices::U, QuantityIndices::V, QuantityIndices::W};
      }
    }();

    // Y+ T+ and Y- T-; the matrices are dense and column major, so [col * Count + row]
    std::array<real, Count> admittedPlus{};
    std::array<real, Count> admittedMinus{};
    for (std::size_t k = 0; k < Count; ++k) {
      for (std::size_t j = 0; j < Count; ++j) {
        admittedPlus[j] += impedanceMatrices.impedance[k * Count + j] *
                           local.faceAlignedValuesPlus[StressIndices[k]];
        admittedMinus[j] += impedanceMatrices.impedanceNeig[k * Count + j] *
                            local.faceAlignedValuesMinus[StressIndices[k]];
      }
    }

    std::array<real, Count> traction{};
    for (std::size_t k = 0; k < Count; ++k) {
      const real rhs = diff(VelocityIndices[k]) + admittedPlus[k] + admittedMinus[k];
      for (std::size_t j = 0; j < Count; ++j) {
        traction[j] += impedanceMatrices.eta[k * Count + j] * rhs;
      }
    }

    local.transientNormalTraction = traction[0];
    local.faceAlignedStress12 = traction[1];
    local.faceAlignedStress13 = traction[2];

    std::array<real, Count> tractionDiff{};
    for (std::size_t k = 0; k < Count; ++k) {
      tractionDiff[k] = traction[k] - local.faceAlignedValuesPlus[StressIndices[k]];
    }

    real normalVelocity = local.faceAlignedValuesPlus[QuantityIndices::U];
    std::array<real, 3> lateralStress{};
    for (std::size_t k = 0; k < Count; ++k) {
      normalVelocity += impedanceMatrices.impedance[k * Count + 0] * tractionDiff[k];
      for (std::size_t j = 0; j < 3; ++j) {
        lateralStress[j] += impedanceMatrices.lateralStress[k * 3 + j] * tractionDiff[k];
      }
    }
    local.faultNormalVelocity = normalVelocity;

    // the stress components which are not part of the fault-normal Riemann problem
    local.faceAlignedStress22 = local.faceAlignedValuesPlus[QuantityIndices::YY] + lateralStress[0];
    local.faceAlignedStress33 = local.faceAlignedValuesPlus[QuantityIndices::ZZ] + lateralStress[1];
    local.faceAlignedStress23 = local.faceAlignedValuesPlus[QuantityIndices::YZ] + lateralStress[2];
  } else {
    const auto& impAndEta = ((local.layer->var<DynamicRupture::ImpAndEta>())[local.ltsId]);
    const real normalDivisor = 1.0 / (impAndEta.zpNeig + impAndEta.zp);
    const real shearDivisor = 1.0 / (impAndEta.zsNeig + impAndEta.zs);

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

    local.faceAlignedStress22 =
        local.faceAlignedValuesPlus[QuantityIndices::YY] + missingSigmaValues;
    local.faceAlignedStress33 =
        local.faceAlignedValuesPlus[QuantityIndices::ZZ] + missingSigmaValues;
    local.faceAlignedStress23 = local.faceAlignedValuesPlus[QuantityIndices::YZ];
  }
}

void ReceiverOutput::updateLocalTractions(LocalInfo& local, real strength, real strengthSlope) {
  const auto component1 = local.iniTraction1 + local.faceAlignedStress12;
  const auto component2 = local.iniTraction2 + local.faceAlignedStress13;
  const auto tracEla = misc::magnitude(component1, component2);

  if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
    // the very solve the friction laws run, so the reconstruction cannot drift away from it: with
    // an anisotropic impedance the slip is not parallel to the trial traction, and the strength
    // follows the fault-normal traction, which follows the slip rate
    const auto& impAndEta = ((local.layer->var<DynamicRupture::ImpAndEta>())[local.ltsId]);
    const auto& impedanceMatrices =
        ((local.layer->var<DynamicRupture::ImpedanceMatrices>())[local.ltsId]);

    const auto solution = friction_law::common::solveSlipRate(
        impAndEta, impedanceMatrices, component1, component2, tracEla, strength, strengthSlope);

    local.slipRateTangent1 = solution.slipRate * solution.direction1;
    local.slipRateTangent2 = solution.slipRate * solution.direction2;

    const auto [tractionUpdate1, tractionUpdate2] = friction_law::common::matmulEta(
        impAndEta, impedanceMatrices, local.slipRateTangent1, local.slipRateTangent2);
    const auto normalUpdate = friction_law::common::matmulEtaNormal(
        impAndEta, impedanceMatrices, local.slipRateTangent1, local.slipRateTangent2);

    local.updatedTraction1 = local.faceAlignedStress12 - tractionUpdate1;
    local.updatedTraction2 = local.faceAlignedStress13 - tractionUpdate2;
    local.transientNormalTraction -= normalUpdate;

    // computeLocalStresses maps the traction of the Riemann problem to the velocity of the Godunov
    // state through the first row of Y+. The friction solve moves that traction, and with an
    // anisotropic admittance the two shear components move the fault-normal velocity as well.
    constexpr std::size_t Count = tensor::Zplus::Shape[0];
    local.faultNormalVelocity -= impedanceMatrices.impedance[0 * Count + 0] * normalUpdate +
                                 impedanceMatrices.impedance[1 * Count + 0] * tractionUpdate1 +
                                 impedanceMatrices.impedance[2 * Count + 0] * tractionUpdate2;
  } else {
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
}

void ReceiverOutput::projectOntoStrikeAndDip(LocalInfo& local,
                                             real alongTangent1,
                                             real alongTangent2,
                                             const std::array<double, 3>& tangent1,
                                             const std::array<double, 3>& tangent2,
                                             const std::array<double, 3>& strike,
                                             const std::array<double, 3>& dip) {
  local.slipRateStrike = static_cast<real>(0.0);
  local.slipRateDip = static_cast<real>(0.0);

  for (size_t i = 0; i < 3; ++i) {
    const real component = alongTangent1 * tangent1[i] + alongTangent2 * tangent2[i];
    local.slipRateStrike += component * strike[i];
    local.slipRateDip += component * dip[i];
  }
}

void ReceiverOutput::computeSlipRate(
    LocalInfo& local,
    [[maybe_unused]] const std::array<real, 6>& rotatedUpdatedStress,
    [[maybe_unused]] const std::array<real, 6>& rotatedStress,
    [[maybe_unused]] const std::array<double, 3>& tangent1,
    [[maybe_unused]] const std::array<double, 3>& tangent2,
    [[maybe_unused]] const std::array<double, 3>& strike,
    [[maybe_unused]] const std::array<double, 3>& dip) {

  if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
    // updateLocalTractions resolves the slip direction along with the magnitude, so all that is
    // left is the rotation onto strike and dip. Recovering the slip rate from the traction
    // difference instead would have to invert the full eta, fault-normal row included, since the
    // shear slip also changes the normal traction.
    projectOntoStrikeAndDip(
        local, local.slipRateTangent1, local.slipRateTangent2, tangent1, tangent2, strike, dip);
  } else {
    // the shear block of eta is a multiple of the identity for every material with an isotropic
    // frame -- poroelasticity included, where the fluid column does not reach the shear rows -- so
    // a scalar is exact and the order of scaling and rotation does not matter
    const auto& impAndEta = ((local.layer->var<DynamicRupture::ImpAndEta>())[local.ltsId]);
    local.slipRateStrike = -impAndEta.invEtaS * (rotatedUpdatedStress[QuantityIndices::XY] -
                                                 rotatedStress[QuantityIndices::XY]);
    local.slipRateDip = -impAndEta.invEtaS * (rotatedUpdatedStress[QuantityIndices::XZ] -
                                              rotatedStress[QuantityIndices::XZ]);
  }
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

real ReceiverOutput::computeRuptureVelocity(const Eigen::Matrix<real, 2, 2>& jacobiT2d,
                                            const LocalInfo& local) {
  const auto* ruptureTime = getCellData<DynamicRupture::RuptureTime>(local);
  real ruptureVelocity = 0.0;

  bool needsUpdate{true};
  for (size_t point = 0; point < misc::NumBoundaryGaussPoints; ++point) {
    if (ruptureTime[point] == 0.0) {
      needsUpdate = false;
    }
  }

  if (needsUpdate) {
    constexpr int NumPoly = ConvergenceOrder - 1;
    constexpr int NumDegFr2d = (NumPoly + 1) * (NumPoly + 2) / 2;
    std::array<double, NumDegFr2d> projectedRT{};
    projectedRT.fill(0.0);

    std::array<double, static_cast<std::size_t>(2 * NumDegFr2d)> phiAtPoint{};
    phiAtPoint.fill(0.0);

    const auto chiTau2dPoints = init::quadpoints::view::create(init::quadpoints::Values);
    const auto weights = init::quadweights::view::create(init::quadweights::Values);

    const auto* rt = getCellData<DynamicRupture::RuptureTime>(local);
    for (size_t jBndGP = 0; jBndGP < misc::NumBoundaryGaussPoints; ++jBndGP) {
      const real chi = seissol::multisim::multisimTranspose(chiTau2dPoints, jBndGP, 0);
      const real tau = seissol::multisim::multisimTranspose(chiTau2dPoints, jBndGP, 1);

      basisFunction::tri_dubiner::evaluatePolynomials(phiAtPoint.data(), chi, tau, NumPoly);

      for (size_t d = 0; d < NumDegFr2d; ++d) {
        projectedRT[d] += weights(jBndGP) * rt[jBndGP] * phiAtPoint[d];
      }
    }
    const auto m2inv = seissol::init::M2inv::view::create(seissol::init::M2inv::Values);
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

std::vector<std::size_t> ReceiverOutput::getOutputVariables() const {
  return {drStorage_->info<DynamicRupture::InitialStressInFaultCS>().index,
          drStorage_->info<DynamicRupture::Mu>().index,
          drStorage_->info<DynamicRupture::RuptureTime>().index,
          drStorage_->info<DynamicRupture::AccumulatedSlipMagnitude>().index,
          drStorage_->info<DynamicRupture::PeakSlipRate>().index,
          drStorage_->info<DynamicRupture::DynStressTime>().index,
          drStorage_->info<DynamicRupture::Slip1>().index,
          drStorage_->info<DynamicRupture::Slip2>().index};
}

} // namespace seissol::dr::output
