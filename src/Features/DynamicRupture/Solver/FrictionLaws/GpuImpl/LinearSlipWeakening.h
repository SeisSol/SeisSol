// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_LINEARSLIPWEAKENING_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_LINEARSLIPWEAKENING_H_

#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/LinearSlipWeakening.h"
#include "Features/DynamicRupture/Solver/Storage.h"

namespace seissol::dr::friction_law::gpu {

template <typename Cfg>
using NoSpecialization = friction_law::NoSpecialization<Cfg>;
template <typename Cfg>
using BiMaterialFault = friction_law::BiMaterialFault<Cfg>;
template <typename Cfg>
using TPApprox = friction_law::TPApprox<Cfg>;

/// The linear slip weakening laws (see friction_law::LinearSlipWeakening) on the device.
template <typename Cfg, typename SpecializationT>
class LinearSlipWeakeningLaw final
    : public BaseFrictionSolver<Cfg, LinearSlipWeakeningLaw<Cfg, SpecializationT>>,
      public LinearSlipWeakening<Cfg, SpecializationT> {
  public:
  explicit LinearSlipWeakeningLaw(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : BaseFrictionSolver<Cfg, LinearSlipWeakeningLaw<Cfg, SpecializationT>>(drParameters) {}

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {
    LinearSlipWeakening<Cfg, SpecializationT>::copySpecificStorageDataToLocal(
        data, layerData, seissol::initializer::AllocationPlace::Device);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_LINEARSLIPWEAKENING_H_
