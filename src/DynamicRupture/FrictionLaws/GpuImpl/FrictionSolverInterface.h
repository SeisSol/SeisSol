// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_

#include "Config.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Typedefs.h"
#include "GeneratedCode/tensor.h"
#include "Memory/Tree/Layer.h"

// A sycl-independent interface is required for interacting with the wp solver
// which, in its turn, is not supposed to know anything about SYCL
namespace seissol::dr::friction_law::gpu {
struct FrictionLawData {
  FrictionLawParameters drParameters;

  const ImpedancesAndEta<Config>* __restrict impAndEta{};
  const ImpedanceMatrices<Config>* __restrict impedanceMatrices{};
  // CS = coordinate system
  const real (*__restrict stressSourceInFaultCS)[6][misc::NumPaddedPoints<Config>]{};
  const real (*__restrict cohesion)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict mu)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict accumulatedSlipMagnitude)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict slip1)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict slip2)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict slipRateMagnitude)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict slipRate1)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict slipRate2)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict ruptureTime)[misc::NumPaddedPoints<Config>]{};
  bool (*__restrict ruptureTimePending)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict peakSlipRate)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict traction1)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict traction2)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict imposedStatePlus)[tensor::QInterpolated<Config>::size()]{};
  real (*__restrict imposedStateMinus)[tensor::QInterpolated<Config>::size()]{};
  DREnergyOutput<Config>* __restrict energyData{};
  const DRGodunovData<Config>* __restrict godunovData{};
  const real (*__restrict stressSourcePressure)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict stressSourceOnset)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict stressSourceRiseTime)[misc::NumPaddedPoints<Config>]{};

  // be careful only for some FLs initialized:
  real (*__restrict dynStressTime)[misc::NumPaddedPoints<Config>]{};
  bool (*__restrict dynStressTimePending)[misc::NumPaddedPoints<Config>]{};

  const real (*__restrict qInterpolatedPlus)[misc::TimeSteps<Config>]
                                            [tensor::QInterpolated<Config>::size()]{};
  const real (*__restrict qInterpolatedMinus)[misc::TimeSteps<Config>]
                                             [tensor::QInterpolated<Config>::size()]{};

  // LSW
  const real (*__restrict dC)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict muS)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict muD)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict forcedRuptureTime)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict regularizedStrength)[misc::NumPaddedPoints<Config>]{};

  // R+S
  const real (*__restrict a)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict sl0)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict stateVariable)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict f0)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict muW)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict b)[misc::NumPaddedPoints<Config>]{};
  bool (*__restrict convergenceInner)[misc::NumPaddedPoints<Config>]{};
  bool (*__restrict convergenceOuter)[misc::NumPaddedPoints<Config>]{};

  // R+S FVW
  const real (*__restrict srW)[misc::NumPaddedPoints<Config>]{};

  // TP
  real (*__restrict temperature)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict pressure)[misc::NumPaddedPoints<Config>]{};
  real (*__restrict theta)[misc::NumTpGridPoints][misc::NumPaddedPoints<Config>]{};
  real (*__restrict sigma)[misc::NumTpGridPoints][misc::NumPaddedPoints<Config>]{};
  const real (*__restrict halfWidthShearZone)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict hydraulicDiffusivity)[misc::NumPaddedPoints<Config>]{};

  // ISR
  const real (*__restrict imposedSlipDirection1)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict imposedSlipDirection2)[misc::NumPaddedPoints<Config>]{};

  // ISR/STF
  const real (*__restrict onsetTime)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict tauS)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict tauR)[misc::NumPaddedPoints<Config>]{};
  const real (*__restrict riseTime)[misc::NumPaddedPoints<Config>]{};
};

class FrictionSolverInterface : public seissol::dr::friction_law::FrictionSolver {
  public:
  explicit FrictionSolverInterface(const FrictionLawParameters& drParameters)
      : seissol::dr::friction_law::FrictionSolver(drParameters) {}
  ~FrictionSolverInterface() override = default;

  seissol::initializer::AllocationPlace allocationPlace() override {
    return seissol::initializer::AllocationPlace::Device;
  }

  static void copyStorageToLocal(FrictionLawData* data, DynamicRupture::Layer& layerData) {
    const seissol::initializer::AllocationPlace place =
        seissol::initializer::AllocationPlace::Device;
    data->impAndEta = layerData.var<DynamicRupture::ImpAndEta>(Config(), place);
    data->impedanceMatrices = layerData.var<DynamicRupture::ImpedanceMatrices>(Config(), place);
    data->stressSourceInFaultCS =
        layerData.var<DynamicRupture::StressSourceInFaultCS>(Config(), place);
    data->mu = layerData.var<DynamicRupture::Mu>(Config(), place);
    data->accumulatedSlipMagnitude =
        layerData.var<DynamicRupture::AccumulatedSlipMagnitude>(Config(), place);
    data->slip1 = layerData.var<DynamicRupture::Slip1>(Config(), place);
    data->slip2 = layerData.var<DynamicRupture::Slip2>(Config(), place);
    data->slipRateMagnitude = layerData.var<DynamicRupture::SlipRateMagnitude>(Config(), place);
    data->slipRate1 = layerData.var<DynamicRupture::SlipRate1>(Config(), place);
    data->slipRate2 = layerData.var<DynamicRupture::SlipRate2>(Config(), place);
    data->ruptureTime = layerData.var<DynamicRupture::RuptureTime>(Config(), place);
    data->ruptureTimePending = layerData.var<DynamicRupture::RuptureTimePending>(Config(), place);
    data->peakSlipRate = layerData.var<DynamicRupture::PeakSlipRate>(Config(), place);
    data->traction1 = layerData.var<DynamicRupture::Traction1>(Config(), place);
    data->traction2 = layerData.var<DynamicRupture::Traction2>(Config(), place);
    data->imposedStatePlus = layerData.var<DynamicRupture::ImposedStatePlus>(Config(), place);
    data->imposedStateMinus = layerData.var<DynamicRupture::ImposedStateMinus>(Config(), place);
    data->energyData = layerData.var<DynamicRupture::DREnergyOutputVar>(Config(), place);
    data->godunovData = layerData.var<DynamicRupture::GodunovData>(Config(), place);
    data->dynStressTime = layerData.var<DynamicRupture::DynStressTime>(Config(), place);
    data->dynStressTimePending =
        layerData.var<DynamicRupture::DynStressTimePending>(Config(), place);
    data->qInterpolatedPlus = layerData.var<DynamicRupture::QInterpolatedPlus>(Config(), place);
    data->qInterpolatedMinus = layerData.var<DynamicRupture::QInterpolatedMinus>(Config(), place);
    data->stressSourcePressure =
        layerData.var<DynamicRupture::StressSourcePressure>(Config(), place);
    data->stressSourceOnset = layerData.var<DynamicRupture::StressSourceOnset>(Config(), place);
    data->stressSourceRiseTime =
        layerData.var<DynamicRupture::StressSourceRiseTime>(Config(), place);
  }

  protected:
  FrictionLawData dataHost_;
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_
