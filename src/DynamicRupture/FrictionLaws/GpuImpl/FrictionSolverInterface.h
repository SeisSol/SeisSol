// SPDX-FileCopyrightText: 2023 SeisSol Group
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

template<typename Cfg>
struct FrictionLawData {
  seissol::initializer::parameters::DRParameters drParameters;

  const ImpedancesAndEta<Real<Cfg>>* __restrict impAndEta{};
  const ImpedanceMatrices<Cfg>* __restrict impedanceMatrices{};
  // CS = coordinate system
  real (*__restrict initialStressInFaultCS)[6][misc::NumPaddedPoints<Cfg>]{};
  const real (*__restrict nucleationStressInFaultCS)[6][misc::NumPaddedPoints<Cfg>]{};
  const real (*__restrict cohesion)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict mu)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict accumulatedSlipMagnitude)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slip1)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slip2)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slipRateMagnitude)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slipRate1)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slipRate2)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict ruptureTime)[misc::NumPaddedPoints<Cfg>]{};
  bool (*__restrict ruptureTimePending)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict peakSlipRate)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict traction1)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict traction2)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict imposedStatePlus)[tensor::QInterpolated<Cfg>::size()]{};
  real (*__restrict imposedStateMinus)[tensor::QInterpolated<Cfg>::size()]{};
  DREnergyOutput<Cfg>* __restrict energyData{};
  const DRGodunovData<Cfg>* __restrict godunovData{};
  real (*__restrict initialPressure)[misc::NumPaddedPoints<Cfg>]{};
  const real (*__restrict nucleationPressure)[misc::NumPaddedPoints<Cfg>]{};

  // be careful only for some FLs initialized:
  real (*__restrict dynStressTime)[misc::NumPaddedPoints<Cfg>]{};
  bool (*__restrict dynStressTimePending)[misc::NumPaddedPoints<Cfg>]{};

  const real (*__restrict qInterpolatedPlus)[Cfg::ConvergenceOrder][tensor::QInterpolated<Cfg>::size()]{};
  const real (*__restrict qInterpolatedMinus)[Cfg::ConvergenceOrder][tensor::QInterpolated<Cfg>::size()]{};

  // LSW
  const real (*__restrict dC)[misc::NumPaddedPoints<Cfg>];
  const real (*__restrict muS)[misc::NumPaddedPoints<Cfg>];
  const real (*__restrict muD)[misc::NumPaddedPoints<Cfg>];
  const real (*__restrict forcedRuptureTime)[misc::NumPaddedPoints<Cfg>];
  real (*__restrict regularizedStrength)[misc::NumPaddedPoints<Cfg>];

  // R+S
  const real (*__restrict a)[misc::NumPaddedPoints<Cfg>];
  const real (*__restrict sl0)[misc::NumPaddedPoints<Cfg>];
  real (*__restrict stateVariable)[misc::NumPaddedPoints<Cfg>];

  // R+S FVW
  const real (*__restrict srW)[misc::NumPaddedPoints<Cfg>];

  // TP
  real (*__restrict temperature)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict pressure)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict theta)[misc::NumTpGridPoints][misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict sigma)[misc::NumTpGridPoints][misc::NumPaddedPoints<Cfg>]{};
  const real (*__restrict halfWidthShearZone)[misc::NumPaddedPoints<Cfg>]{};
  const real (*__restrict hydraulicDiffusivity)[misc::NumPaddedPoints<Cfg>]{};

  // ISR
  const real (*__restrict imposedSlipDirection1)[misc::NumPaddedPoints<Cfg>];
  const real (*__restrict imposedSlipDirection2)[misc::NumPaddedPoints<Cfg>];

  // ISR/STF
  const real (*__restrict onsetTime)[misc::NumPaddedPoints<Cfg>];
  const real (*__restrict tauS)[misc::NumPaddedPoints<Cfg>];
  const real (*__restrict tauR)[misc::NumPaddedPoints<Cfg>];
  const real (*__restrict riseTime)[misc::NumPaddedPoints<Cfg>];
};

class FrictionSolverInterface : public seissol::dr::friction_law::FrictionSolver {
  public:
  explicit FrictionSolverInterface(seissol::initializer::parameters::DRParameters* drParameters)
      : seissol::dr::friction_law::FrictionSolver(drParameters) {}
  ~FrictionSolverInterface() override = default;

  seissol::initializer::AllocationPlace allocationPlace() override {
    return seissol::initializer::AllocationPlace::Device;
  }

  static void copyStorageToLocal(FrictionLawData* data, DynamicRupture::Layer& layerData) {
    const seissol::initializer::AllocationPlace place =
        seissol::initializer::AllocationPlace::Device;
    data->impAndEta = layerData.var<DynamicRupture::ImpAndEta>(place);
    data->impedanceMatrices = layerData.var<DynamicRupture::ImpedanceMatrices>(place);
    data->initialStressInFaultCS = layerData.var<DynamicRupture::InitialStressInFaultCS>(place);
    data->nucleationStressInFaultCS =
        layerData.var<DynamicRupture::NucleationStressInFaultCS>(place);
    data->mu = layerData.var<DynamicRupture::Mu>(place);
    data->accumulatedSlipMagnitude = layerData.var<DynamicRupture::AccumulatedSlipMagnitude>(place);
    data->slip1 = layerData.var<DynamicRupture::Slip1>(place);
    data->slip2 = layerData.var<DynamicRupture::Slip2>(place);
    data->slipRateMagnitude = layerData.var<DynamicRupture::SlipRateMagnitude>(place);
    data->slipRate1 = layerData.var<DynamicRupture::SlipRate1>(place);
    data->slipRate2 = layerData.var<DynamicRupture::SlipRate2>(place);
    data->ruptureTime = layerData.var<DynamicRupture::RuptureTime>(place);
    data->ruptureTimePending = layerData.var<DynamicRupture::RuptureTimePending>(place);
    data->peakSlipRate = layerData.var<DynamicRupture::PeakSlipRate>(place);
    data->traction1 = layerData.var<DynamicRupture::Traction1>(place);
    data->traction2 = layerData.var<DynamicRupture::Traction2>(place);
    data->imposedStatePlus = layerData.var<DynamicRupture::ImposedStatePlus>(place);
    data->imposedStateMinus = layerData.var<DynamicRupture::ImposedStateMinus>(place);
    data->energyData = layerData.var<DynamicRupture::DREnergyOutputVar>(place);
    data->godunovData = layerData.var<DynamicRupture::GodunovData>(place);
    data->dynStressTime = layerData.var<DynamicRupture::DynStressTime>(place);
    data->dynStressTimePending = layerData.var<DynamicRupture::DynStressTimePending>(place);
    data->qInterpolatedPlus = layerData.var<DynamicRupture::QInterpolatedPlus>(place);
    data->qInterpolatedMinus = layerData.var<DynamicRupture::QInterpolatedMinus>(place);
    data->initialPressure = layerData.var<DynamicRupture::InitialPressure>(place);
    data->nucleationPressure = layerData.var<DynamicRupture::NucleationPressure>(place);
  }

  protected:
  FrictionLawData dataHost;
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_
