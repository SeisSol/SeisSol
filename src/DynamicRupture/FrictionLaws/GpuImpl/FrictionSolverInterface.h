// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "Initializer/Parameters/DRParameters.h"
#include <Memory/Tree/Layer.h>

// A sycl-independent interface is required for interacting with the wp solver
// which, in its turn, is not supposed to know anything about SYCL
namespace seissol::dr::friction_law::gpu {
struct FrictionLawData {
  real deltaT[ConvergenceOrder] = {};
  real sumDt{};

  seissol::initializer::parameters::DRParameters drParameters;

  ImpedancesAndEta* impAndEta{};
  ImpedanceMatrices* impedanceMatrices{};
  real mFullUpdateTime{};
  // CS = coordinate system
  real (*initialStressInFaultCS)[misc::NumPaddedPoints][6]{};
  real (*nucleationStressInFaultCS)[misc::NumPaddedPoints][6]{};
  real (*cohesion)[misc::NumPaddedPoints]{};
  real (*mu)[misc::NumPaddedPoints]{};
  real (*accumulatedSlipMagnitude)[misc::NumPaddedPoints]{};
  real (*slip1)[misc::NumPaddedPoints]{};
  real (*slip2)[misc::NumPaddedPoints]{};
  real (*slipRateMagnitude)[misc::NumPaddedPoints]{};
  real (*slipRate1)[misc::NumPaddedPoints]{};
  real (*slipRate2)[misc::NumPaddedPoints]{};
  real (*ruptureTime)[misc::NumPaddedPoints]{};
  bool (*ruptureTimePending)[misc::NumPaddedPoints]{};
  real (*peakSlipRate)[misc::NumPaddedPoints]{};
  real (*traction1)[misc::NumPaddedPoints]{};
  real (*traction2)[misc::NumPaddedPoints]{};
  real (*imposedStatePlus)[tensor::QInterpolated::size()]{};
  real (*imposedStateMinus)[tensor::QInterpolated::size()]{};
  real spaceWeights[misc::NumPaddedPoints]{};
  DREnergyOutput* energyData{};
  DRGodunovData* godunovData{};
  real (*initialPressure)[misc::NumPaddedPoints]{};
  real (*nucleationPressure)[misc::NumPaddedPoints]{};

  // be careful only for some FLs initialized:
  real (*dynStressTime)[misc::NumPaddedPoints]{};
  bool (*dynStressTimePending)[misc::NumPaddedPoints]{};

  real (*qInterpolatedPlus)[ConvergenceOrder][tensor::QInterpolated::size()]{};
  real (*qInterpolatedMinus)[ConvergenceOrder][tensor::QInterpolated::size()]{};

  // LSW
  real (*dC)[misc::NumPaddedPoints];
  real (*muS)[misc::NumPaddedPoints];
  real (*muD)[misc::NumPaddedPoints];
  real (*forcedRuptureTime)[misc::NumPaddedPoints];
  real (*regularisedStrength)[misc::NumPaddedPoints];

  // R+S
  real (*a)[misc::NumPaddedPoints];
  real (*sl0)[misc::NumPaddedPoints];
  real (*stateVariable)[misc::NumPaddedPoints];

  // R+S FVW
  real (*srW)[misc::NumPaddedPoints];
};

class FrictionSolverInterface : public seissol::dr::friction_law::FrictionSolver {
  public:
  explicit FrictionSolverInterface(seissol::initializer::parameters::DRParameters* drParameters)
      : seissol::dr::friction_law::FrictionSolver(drParameters) {}
  ~FrictionSolverInterface() override = default;

  virtual void initSyclQueue() = 0;
  void setMaxClusterSize(size_t size) { maxClusterSize = size; }
  virtual void allocateAuxiliaryMemory() = 0;
  virtual void copyStaticDataToDevice() = 0;

  seissol::initializer::AllocationPlace allocationPlace() override {
    return seissol::initializer::AllocationPlace::Device;
  }

  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* dynRup,
                                 real fullUpdateTime) {
    const seissol::initializer::AllocationPlace place =
        seissol::initializer::AllocationPlace::Device;
    data->impAndEta = layerData.var(dynRup->impAndEta, place);
    data->impedanceMatrices = layerData.var(dynRup->impedanceMatrices, place);
    data->initialStressInFaultCS = layerData.var(dynRup->initialStressInFaultCS, place);
    data->nucleationStressInFaultCS = layerData.var(dynRup->nucleationStressInFaultCS, place);
    data->mu = layerData.var(dynRup->mu, place);
    data->accumulatedSlipMagnitude = layerData.var(dynRup->accumulatedSlipMagnitude, place);
    data->slip1 = layerData.var(dynRup->slip1, place);
    data->slip2 = layerData.var(dynRup->slip2, place);
    data->slipRateMagnitude = layerData.var(dynRup->slipRateMagnitude, place);
    data->slipRate1 = layerData.var(dynRup->slipRate1, place);
    data->slipRate2 = layerData.var(dynRup->slipRate2, place);
    data->ruptureTime = layerData.var(dynRup->ruptureTime, place);
    data->ruptureTimePending = layerData.var(dynRup->ruptureTimePending, place);
    data->peakSlipRate = layerData.var(dynRup->peakSlipRate, place);
    data->traction1 = layerData.var(dynRup->traction1, place);
    data->traction2 = layerData.var(dynRup->traction2, place);
    data->imposedStatePlus = layerData.var(dynRup->imposedStatePlus, place);
    data->imposedStateMinus = layerData.var(dynRup->imposedStateMinus, place);
    data->energyData = layerData.var(dynRup->drEnergyOutput, place);
    data->godunovData = layerData.var(dynRup->godunovData, place);
    data->mFullUpdateTime = fullUpdateTime;
    data->dynStressTime = layerData.var(dynRup->dynStressTime, place);
    data->dynStressTimePending = layerData.var(dynRup->dynStressTimePending, place);
    data->qInterpolatedPlus = layerData.var(dynRup->qInterpolatedPlus, place);
    data->qInterpolatedMinus = layerData.var(dynRup->qInterpolatedMinus, place);
    data->initialPressure = layerData.var(dynRup->initialPressure, place);
    data->nucleationPressure = layerData.var(dynRup->nucleationPressure, place);
  }

  protected:
  size_t maxClusterSize{};
  FrictionLawData dataHost;
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_
