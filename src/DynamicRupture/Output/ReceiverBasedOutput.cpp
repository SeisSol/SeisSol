// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ReceiverBasedOutput.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Output/ImposedSlipRates.h"
#include "DynamicRupture/Output/LinearSlipWeakening.h"
#include "DynamicRupture/Output/LinearSlipWeakeningBimaterial.h"
#include "DynamicRupture/Output/NoFault.h"
#include "DynamicRupture/Output/RateAndState.h"
#include "DynamicRupture/Output/RateAndStateThermalPressurization.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Common.h"
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
#include <cstdint>
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

std::vector<std::size_t> ReceiverOutput::getOutputVariables() const {
  return {drStorage_->info<DynamicRupture::StressSourceInFaultCS>().index,
          drStorage_->info<DynamicRupture::Mu>().index,
          drStorage_->info<DynamicRupture::RuptureTime>().index,
          drStorage_->info<DynamicRupture::AccumulatedSlipMagnitude>().index,
          drStorage_->info<DynamicRupture::PeakSlipRate>().index,
          drStorage_->info<DynamicRupture::DynStressTime>().index,
          drStorage_->info<DynamicRupture::Slip1>().index,
          drStorage_->info<DynamicRupture::Slip2>().index};
}

template <typename Derived>
template <typename Cfg>
void ReceiverOutputImpl<Derived>::getDofs(const Real<Cfg>*(&derivatives), std::size_t meshId) {
  const auto position = wpBackmap_->get(meshId);
  auto& layer = wpStorage_->layer(position.color);
  // get DOFs from 0th derivatives
  assert(
      layer.var<LTS::CellInformation>()[position.cell].ltsSetup.hasBuffer(BufferType::Derivatives));
  // the cells next to a fault face compute in the configuration of the face
  assert(layer.getIdentifier().config == configIdOf<Cfg>());

  derivatives = layer.var<LTS::Derivatives>(Cfg())[position.cell];
}

template <typename Derived>
template <typename Cfg>
void ReceiverOutputImpl<Derived>::getNeighborDofs(const Real<Cfg>*(&derivatives),
                                                  std::size_t meshId,
                                                  std::size_t side) {
  const auto position = wpBackmap_->get(meshId);
  auto& layer = wpStorage_->layer(position.color);

  derivatives = static_cast<const Real<Cfg>*>(layer.var<LTS::FaceNeighbors>()[position.cell][side]);
  assert(derivatives != nullptr);
}

template <typename Derived>
void ReceiverOutputImpl<Derived>::calcFaultOutput(
    seissol::initializer::parameters::OutputType outputType,
    seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
    const std::shared_ptr<ReceiverOutputData>& outputData,
    parallel::runtime::StreamRuntime& runtime,
    double stateTime,
    double time,
    double dt,
    double indt) {
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

  forEachConfig([&](auto cfg) {
    using Cfg = decltype(cfg);
    this->template calcFaultOutputOfConfig<Cfg>(
        outputType, slipRateOutputType, outputData, runtime, stateTime, dt, indt);
  });

  if (outputType == seissol::initializer::parameters::OutputType::AtPickpoint) {
    outputData->cachedTime[outputData->currentCacheLevel] = time;
    outputData->currentCacheLevel += 1;
  }
}

template <typename Derived>
template <typename Cfg>
void ReceiverOutputImpl<Derived>::calcFaultOutputOfConfig(
    seissol::initializer::parameters::OutputType outputType,
    seissol::initializer::parameters::SlipRateOutputType slipRateOutputType,
    const std::shared_ptr<ReceiverOutputData>& outputData,
    parallel::runtime::StreamRuntime& runtime,
    double stateTime,
    double dt,
    double indt) {
  using real = Real<Cfg>;

  const size_t level = (outputType == seissol::initializer::parameters::OutputType::AtPickpoint)
                           ? outputData->currentCacheLevel
                           : 0;
  const auto& faultInfos = meshReader_->getFault();

  // the friction solve advances in the sub intervals of the time quadrature; the stored friction
  // state belongs to the last of them
  const auto frictionTime = seissol::dr::friction_law::FrictionSolver::computeDeltaT<Cfg>(
      seissol::quadrature::ShiftedGaussLegendre(Cfg::ConvergenceOrder, 0, dt).first);

  const auto timeCoeffs = kernels::timeBasis<Cfg>().point(indt, dt);

  auto& callRuntime =
      outputData->extraRuntime.has_value() ? outputData->extraRuntime.value() : runtime;

  const auto handler = [this,
                        outputData,
                        &faultInfos,
                        outputType,
                        slipRateOutputType,
                        level,
                        timeCoeffs,
                        stateTime,
                        frictionTime](std::size_t faceId) {
    constexpr auto Variant = configIdOf<Cfg>();

    const auto& topology = outputData->topology;
    const auto& outFace = topology.faces[faceId];
    auto& faceLayer = drStorage_->layer(outFace.position.color);
    if (faceLayer.getIdentifier().config != Variant) {
      // the face is output by the pass of its own configuration
      return;
    }
    const auto& transform = std::get<FaceTransform<Cfg>>(outFace.transform);

    alignas(Alignment) real dofsPlus[tensor::Q<Cfg>::size()]{};
    alignas(Alignment) real dofsMinus[tensor::Q<Cfg>::size()]{};

    alignas(Alignment) real faceAlignedValuesPlus[tensor::QAtPoint<Cfg>::size()]{};
    alignas(Alignment) real faceAlignedValuesMinus[tensor::QAtPoint<Cfg>::size()]{};

    const auto faceIndex = outFace.faultFaceIndex;

    LocalInfo<Cfg> local{};

    local.layer = &faceLayer;
    local.ltsId = outFace.position.cell;
    local.faceId = faceId;
    local.state = outputData.get();
    local.time = stateTime;
    local.deltaT = frictionTime.deltaT.back();
    local.printWarning = &this->printRSFWarning_;

    local.waveSpeedsPlus =
        &((local.layer->template var<DynamicRupture::WaveSpeedsPlus>())[local.ltsId]);
    local.waveSpeedsMinus =
        &((local.layer->template var<DynamicRupture::WaveSpeedsMinus>())[local.ltsId]);
    const auto& faultInfo = faultInfos[faceIndex];

    if (outputType == initializer::parameters::OutputType::Elementwise) {
      std::memcpy(dofsPlus,
                  local.layer->template var<DynamicRupture::TimeDofsPlus>(Cfg())[local.ltsId],
                  sizeof(dofsPlus));
      std::memcpy(dofsMinus,
                  local.layer->template var<DynamicRupture::TimeDofsMinus>(Cfg())[local.ltsId],
                  sizeof(dofsMinus));
    } else {
      // only interpolate for the on-fault receivers
      const real* stePlus = nullptr;
      const real* steMinus = nullptr;

      if constexpr (isDeviceOn()) {
        stePlus =
            static_cast<const real*>(outputData->deviceDataCollector->get(outFace.deviceDataPlus));
        steMinus =
            static_cast<const real*>(outputData->deviceDataCollector->get(outFace.deviceDataMinus));
      } else {
        getDofs<Cfg>(stePlus, faultInfo.element.value());
        if (faultInfo.neighborElement.hasValue()) {
          getDofs<Cfg>(steMinus, faultInfo.neighborElement.value());
        } else {
          getNeighborDofs<Cfg>(steMinus, faultInfo.element.value(), faultInfo.side);
        }
      }

      kernels::Time<Cfg> timeKernel;
      timeKernel.evaluate(timeCoeffs.data(), stePlus, dofsPlus);
      timeKernel.evaluate(timeCoeffs.data(), steMinus, dofsMinus);
    }

    // the rotations and the interpolation frame are properties of the face, so both kernels are
    // configured once here and only fed with per-point basis functions below
    const auto& normal = outFace.faultDirections.faceNormal;
    const auto& tangent1 = outFace.faultDirections.tangent1;
    const auto& tangent2 = outFace.faultDirections.tangent2;
    const auto& strike = outFace.faultDirections.strike;
    const auto& dip = outFace.faultDirections.dip;
    const auto& jacobiT2d = transform.jacobianT2d;

    const auto sourceCount = stressSourceCount(*drParameters_);
    const auto* stressSources =
        local.layer->template var<DynamicRupture::StressSourceInFaultCS>(Cfg());
    const auto* stressSourceOnset =
        local.layer->template var<DynamicRupture::StressSourceOnset>(Cfg());
    const auto* stressSourceRiseTime =
        local.layer->template var<DynamicRupture::StressSourceRiseTime>(Cfg());

    runtime::dynamicRupture::kernel::evaluateFaceAlignedDOFSAtPoint kernel;
    kernel.Tinv = runtime::init::Tinv::view(Variant, transform.glbToFaceAlignedData.data());

    runtime::dynamicRupture::kernel::rotateInitStress alignAlongDipAndStrikeKernel;
    alignAlongDipAndStrikeKernel.stressRotationMatrix = runtime::init::stressRotationMatrix::view(
        Variant, transform.stressGlbToDipStrikeAligned.data());
    alignAlongDipAndStrikeKernel.reducedFaceAlignedMatrix =
        runtime::init::reducedFaceAlignedMatrix::view(Variant,
                                                      transform.stressFaceAlignedToGlb.data());

    for (const auto pointId : topology.pointsOf(faceId)) {
      const auto& outPoint = topology.points[pointId];
      const auto& basisFunctions = std::get<PlusMinusBasisFunctions<Cfg>>(outPoint.basisFunctions);

      kernel.Q = runtime::init::Q::view(Variant, dofsPlus);
      kernel.basisFunctionsAtPoint =
          runtime::init::basisFunctionsAtPoint::view(Variant, basisFunctions.plusSide.data());
      kernel.QAtPoint = runtime::init::QAtPoint::view(Variant, faceAlignedValuesPlus);
      kernel.execute(Variant);

      kernel.Q = runtime::init::Q::view(Variant, dofsMinus);
      kernel.basisFunctionsAtPoint =
          runtime::init::basisFunctionsAtPoint::view(Variant, basisFunctions.minusSide.data());
      kernel.QAtPoint = runtime::init::QAtPoint::view(Variant, faceAlignedValuesMinus);
      kernel.execute(Variant);

      local.nearestGpIndex = static_cast<int>(outPoint.nearestGpIndex);
      local.nearestInternalGpIndex = static_cast<int>(outPoint.nearestInternalGpIndex);

      for (const auto i : topology.receiversOf(pointId)) {
        local.index = i;
        local.fusedIndex = outputData->receivers[i].simIndex;

        assert(outputData->receivers[i].isInside == true &&
               "a receiver is not within any tetrahedron adjacent to a fault");

        local.gpIndex = outputData->receivers[i].gpIndex;
        local.internalGpIndexFused = outputData->receivers[i].internalGpIndexFused;

        local.frictionCoefficient = getCellData<DynamicRupture::Mu>(local)[local.gpIndex];
        local.stateVariable = derived().computeStateVariable(local);

        // the whole tensor, since the total traction output rotates it
        const auto initialStress =
            stressAtTime<Cfg>(&stressSources[local.ltsId * sourceCount],
                              &stressSourceRiseTime[local.ltsId * sourceCount],
                              &stressSourceOnset[local.ltsId * sourceCount],
                              sourceCount,
                              static_cast<std::uint32_t>(local.gpIndex),
                              static_cast<real>(local.time));

        local.iniTraction1 = initialStress[QuantityIndices::XY];
        local.iniTraction2 = initialStress[QuantityIndices::XZ];
        local.iniNormalTraction = initialStress[QuantityIndices::XX];
        local.fluidPressure = derived().computeFluidPressure(local);

        for (size_t j = 0; j < tensor::QAtPoint<Cfg>::Shape[seissol::multisim::BasisDim<Cfg>];
             ++j) {
          local.faceAlignedValuesPlus[j] =
              faceAlignedValuesPlus[j * Cfg::NumSimulations + local.fusedIndex];
          local.faceAlignedValuesMinus[j] =
              faceAlignedValuesMinus[j * Cfg::NumSimulations + local.fusedIndex];
        }

        derived().handleNonConvergence(local);

        this->computeLocalStresses(local);
        const real strength = derived().computeLocalStrength(local);
        const real strengthSlope = derived().computeLocalStrengthSlope(local);
        updateLocalTractions(local, strength, strengthSlope);

        std::array<real, 6> updatedStress{};
        updatedStress[QuantityIndices::XX] = local.transientNormalTraction;
        updatedStress[QuantityIndices::YY] = local.faceAlignedStress22;
        updatedStress[QuantityIndices::ZZ] = local.faceAlignedStress33;
        updatedStress[QuantityIndices::XY] = local.updatedTraction1;
        updatedStress[QuantityIndices::YZ] = local.faceAlignedStress23;
        updatedStress[QuantityIndices::XZ] = local.updatedTraction2;

        alignAlongDipAndStrikeKernel.initialStress =
            runtime::init::initialStress::view(Variant, updatedStress.data());
        std::array<real, 6> rotatedUpdatedStress{};
        alignAlongDipAndStrikeKernel.rotatedStress =
            runtime::init::rotatedStress::view(Variant, rotatedUpdatedStress.data());
        alignAlongDipAndStrikeKernel.execute(Variant);

        std::array<real, 6> stress{};
        stress[QuantityIndices::XX] = local.transientNormalTraction;
        stress[QuantityIndices::YY] = local.faceAlignedStress22;
        stress[QuantityIndices::ZZ] = local.faceAlignedStress33;
        stress[QuantityIndices::XY] = local.faceAlignedStress12;
        stress[QuantityIndices::YZ] = local.faceAlignedStress23;
        stress[QuantityIndices::XZ] = local.faceAlignedStress13;

        alignAlongDipAndStrikeKernel.initialStress =
            runtime::init::initialStress::view(Variant, stress.data());
        std::array<real, 6> rotatedStress{};
        alignAlongDipAndStrikeKernel.rotatedStress =
            runtime::init::rotatedStress::view(Variant, rotatedStress.data());
        alignAlongDipAndStrikeKernel.execute(Variant);

        switch (slipRateOutputType) {
        case seissol::initializer::parameters::SlipRateOutputType::TractionsAndFailure: {
          derived().computeSlipRate(
              local, rotatedUpdatedStress, rotatedStress, tangent1, tangent2, strike, dip);
          break;
        }
        case seissol::initializer::parameters::SlipRateOutputType::VelocityDifference: {
          computeSlipRate(local, tangent1, tangent2, strike, dip);
          break;
        }
        }

        derived().template adjustRotatedUpdatedStress<Cfg>(rotatedUpdatedStress, rotatedStress);

        auto& slipRate = std::get<VariableID::SlipRate>(outputData->vars);
        if (slipRate.isActive) {
          slipRate(DirectionID::Strike, level, i) = local.slipRateStrike;
          slipRate(DirectionID::Dip, level, i) = local.slipRateDip;
        }

        auto& transientTractions = std::get<VariableID::TransientTractions>(outputData->vars);
        if (transientTractions.isActive) {
          transientTractions(DirectionID::Strike, level, i) =
              rotatedUpdatedStress[QuantityIndices::XY];
          transientTractions(DirectionID::Dip, level, i) =
              rotatedUpdatedStress[QuantityIndices::XZ];
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
          std::array<real, tensor::initialStress<Cfg>::size()> unrotatedInitStress{};
          std::array<real, tensor::rotatedStress<Cfg>::size()> rotatedInitStress{};
          for (std::size_t stressVar = 0; stressVar < unrotatedInitStress.size(); ++stressVar) {
            unrotatedInitStress[stressVar] = initialStress[stressVar];
          }
          alignAlongDipAndStrikeKernel.initialStress =
              runtime::init::initialStress::view(Variant, unrotatedInitStress.data());
          alignAlongDipAndStrikeKernel.rotatedStress =
              runtime::init::rotatedStress::view(Variant, rotatedInitStress.data());
          alignAlongDipAndStrikeKernel.execute(Variant);

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
          CoordinateT crossProduct = {0.0, 0.0, 0.0};
          MeshTools::cross(strike, tangent1, crossProduct);

          const double cos1t = MeshTools::dot(strike, tangent1);
          const double scalarProd = MeshTools::dot(crossProduct, normal);

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
        derived().outputSpecifics(outputData, local, level, i);
      }
    }
  };

  callRuntime.enqueueLoop(outputData->topology.faceCount(), handler);
}

template <typename Derived>
template <typename Cfg>
void ReceiverOutputImpl<Derived>::computeLocalStresses(LocalInfo<Cfg>& local) {
  using real = Real<Cfg>;

  auto diff = [&local](int i) {
    return local.faceAlignedValuesMinus[i] - local.faceAlignedValuesPlus[i];
  };

  // named positively: the matrices below are only filled for the materials that go through the
  // general branch of initializeDynamicRuptureMatrices, and reading them for anything else would
  // reconstruct the Godunov state from zeros
  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic ||
                model::MaterialOf<Cfg>::Type == model::MaterialType::Poroelastic) {
    // Anisotropy couples the fault-normal and the two tangential directions, poroelasticity adds
    // the fluid pressure as a fourth interface variable. In both cases the Godunov state has to be
    // reconstructed with the full matrix -- exactly as
    // common::precomputeStressFromQInterpolated does it for the solver:
    //
    //   T*   = eta * (Y+ T+ + Y- T- + (v- - v+))
    //   v*_+ = v+ + Y+ (T* - T+)
    const auto& impedanceMatrices =
        ((local.layer->template var<DynamicRupture::ImpedanceMatrices>(Cfg()))[local.ltsId]);

    constexpr std::size_t Count =
        model::MaterialOf<Cfg>::Type == model::MaterialType::Poroelastic ? 4 : 3;
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
    const auto& impAndEta =
        ((local.layer->template var<DynamicRupture::ImpAndEta>(Cfg()))[local.ltsId]);
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

template <typename Derived>
template <typename Cfg>
void ReceiverOutputImpl<Derived>::updateLocalTractions(LocalInfo<Cfg>& local,
                                                       Real<Cfg> strength,
                                                       Real<Cfg> strengthSlope) {
  const auto component1 = local.iniTraction1 + local.faceAlignedStress12;
  const auto component2 = local.iniTraction2 + local.faceAlignedStress13;
  const auto tracEla = misc::magnitude(component1, component2);

  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
    // the very solve the friction laws run, so the reconstruction cannot drift away from it: with
    // an anisotropic impedance the slip is not parallel to the trial traction, and the strength
    // follows the fault-normal traction, which follows the slip rate
    const auto& impAndEta =
        ((local.layer->template var<DynamicRupture::ImpAndEta>(Cfg()))[local.ltsId]);
    const auto& impedanceMatrices =
        ((local.layer->template var<DynamicRupture::ImpedanceMatrices>(Cfg()))[local.ltsId]);

    const auto solution = friction_law::common::solveSlipRate<Cfg>(
        impAndEta, impedanceMatrices, component1, component2, tracEla, strength, strengthSlope);

    local.slipRateTangent1 = solution.slipRate * solution.direction1;
    local.slipRateTangent2 = solution.slipRate * solution.direction2;

    const auto [tractionUpdate1, tractionUpdate2] = friction_law::common::matmulEta<Cfg>(
        impAndEta, impedanceMatrices, local.slipRateTangent1, local.slipRateTangent2);
    const auto normalUpdate = friction_law::common::matmulEtaNormal<Cfg>(
        impAndEta, impedanceMatrices, local.slipRateTangent1, local.slipRateTangent2);

    local.updatedTraction1 = local.faceAlignedStress12 - tractionUpdate1;
    local.updatedTraction2 = local.faceAlignedStress13 - tractionUpdate2;
    local.transientNormalTraction -= normalUpdate;

    // computeLocalStresses maps the traction of the Riemann problem to the velocity of the Godunov
    // state through the first row of Y+. The friction solve moves that traction, and with an
    // anisotropic admittance the two shear components move the fault-normal velocity as well.
    constexpr std::size_t Count = tensor::Zplus<Cfg>::Shape[0];
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

template <typename Derived>
template <typename Cfg>
void ReceiverOutputImpl<Derived>::projectOntoStrikeAndDip(LocalInfo<Cfg>& local,
                                                          Real<Cfg> alongTangent1,
                                                          Real<Cfg> alongTangent2,
                                                          const std::array<double, 3>& tangent1,
                                                          const std::array<double, 3>& tangent2,
                                                          const std::array<double, 3>& strike,
                                                          const std::array<double, 3>& dip) {
  local.slipRateStrike = static_cast<Real<Cfg>>(0.0);
  local.slipRateDip = static_cast<Real<Cfg>>(0.0);

  for (size_t i = 0; i < 3; ++i) {
    const Real<Cfg> component = alongTangent1 * tangent1[i] + alongTangent2 * tangent2[i];
    local.slipRateStrike += component * strike[i];
    local.slipRateDip += component * dip[i];
  }
}

template <typename Derived>
template <typename Cfg>
void ReceiverOutputImpl<Derived>::computeSlipRate(
    LocalInfo<Cfg>& local,
    [[maybe_unused]] const std::array<Real<Cfg>, 6>& rotatedUpdatedStress,
    [[maybe_unused]] const std::array<Real<Cfg>, 6>& rotatedStress,
    [[maybe_unused]] const std::array<double, 3>& tangent1,
    [[maybe_unused]] const std::array<double, 3>& tangent2,
    [[maybe_unused]] const std::array<double, 3>& strike,
    [[maybe_unused]] const std::array<double, 3>& dip) {

  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
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
    const auto& impAndEta =
        ((local.layer->template var<DynamicRupture::ImpAndEta>(Cfg()))[local.ltsId]);
    local.slipRateStrike = -impAndEta.invEtaS * (rotatedUpdatedStress[QuantityIndices::XY] -
                                                 rotatedStress[QuantityIndices::XY]);
    local.slipRateDip = -impAndEta.invEtaS * (rotatedUpdatedStress[QuantityIndices::XZ] -
                                              rotatedStress[QuantityIndices::XZ]);
  }
}

template <typename Derived>
template <typename Cfg>
void ReceiverOutputImpl<Derived>::computeSlipRate(LocalInfo<Cfg>& local,
                                                  const std::array<double, 3>& tangent1,
                                                  const std::array<double, 3>& tangent2,
                                                  const std::array<double, 3>& strike,
                                                  const std::array<double, 3>& dip) {
  local.slipRateStrike = static_cast<Real<Cfg>>(0.0);
  local.slipRateDip = static_cast<Real<Cfg>>(0.0);

  for (size_t i = 0; i < 3; ++i) {
    const Real<Cfg> factorMinus = (local.faceAlignedValuesMinus[QuantityIndices::V] * tangent1[i] +
                                   local.faceAlignedValuesMinus[QuantityIndices::W] * tangent2[i]);

    const Real<Cfg> factorPlus = (local.faceAlignedValuesPlus[QuantityIndices::V] * tangent1[i] +
                                  local.faceAlignedValuesPlus[QuantityIndices::W] * tangent2[i]);

    local.slipRateStrike += (factorMinus - factorPlus) * strike[i];
    local.slipRateDip += (factorMinus - factorPlus) * dip[i];
  }
}

template <typename Derived>
template <typename Cfg>
Real<Cfg> ReceiverOutputImpl<Derived>::computeRuptureVelocity(
    const Eigen::Matrix<Real<Cfg>, 2, 2>& jacobiT2d, const LocalInfo<Cfg>& local) {
  using real = Real<Cfg>;
  const auto* ruptureTime = getCellData<DynamicRupture::RuptureTime>(local);
  real ruptureVelocity = 0.0;

  bool needsUpdate{true};
  for (size_t point = 0; point < misc::NumBoundaryGaussPoints<Cfg>; ++point) {
    if (ruptureTime[point * Cfg::NumSimulations + local.fusedIndex] == 0.0) {
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

    const auto chiTau2dPoints = init::quadpoints<Cfg>::view::create(init::quadpoints<Cfg>::Values);
    const auto weights = init::quadweights<Cfg>::view::create(init::quadweights<Cfg>::Values);

    const auto* rt = getCellData<DynamicRupture::RuptureTime>(local);
    for (size_t jBndGP = 0; jBndGP < misc::NumBoundaryGaussPoints<Cfg>; ++jBndGP) {
      const real chi = seissol::multisim::multisimTranspose<Cfg>(chiTau2dPoints, jBndGP, 0);
      const real tau = seissol::multisim::multisimTranspose<Cfg>(chiTau2dPoints, jBndGP, 1);

      basisFunction::tri_dubiner::evaluatePolynomials(phiAtPoint.data(), chi, tau, NumPoly);

      for (size_t d = 0; d < NumDegFr2d; ++d) {
        projectedRT[d] +=
            weights(jBndGP) * rt[jBndGP * Cfg::NumSimulations + local.fusedIndex] * phiAtPoint[d];
      }
    }
    const auto m2inv = seissol::init::M2inv<Cfg>::view::create(seissol::init::M2inv<Cfg>::Values);
    for (size_t d = 0; d < NumDegFr2d; ++d) {
      projectedRT[d] *= m2inv(d, d);
    }

    const real chi =
        seissol::multisim::multisimTranspose<Cfg>(chiTau2dPoints, local.nearestInternalGpIndex, 0);
    const real tau =
        seissol::multisim::multisimTranspose<Cfg>(chiTau2dPoints, local.nearestInternalGpIndex, 1);
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

template class ReceiverOutputImpl<ImposedSlipRates>;
template class ReceiverOutputImpl<LinearSlipWeakening>;
template class ReceiverOutputImpl<LinearSlipWeakeningBimaterial>;
template class ReceiverOutputImpl<NoFault>;
template class ReceiverOutputImpl<RateAndState>;
template class ReceiverOutputImpl<RateAndStateThermalPressurization>;

} // namespace seissol::dr::output
