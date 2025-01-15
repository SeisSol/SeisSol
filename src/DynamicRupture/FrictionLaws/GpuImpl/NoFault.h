// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"

namespace seissol::dr::friction_law::gpu {

class NoFault : public BaseFrictionSolver<NoFault> {
  public:
  NoFault(seissol::initializer::parameters::DRParameters* drParameters)
      : BaseFrictionSolver<NoFault>(drParameters) {};

  void copySpecificLtsDataTreeToLocal(seissol::initializer::Layer& layerData,
                                      const seissol::initializer::DynamicRupture* const dynRup,
                                      real fullUpdateTime) override {}

  void updateFrictionAndSlip(unsigned timeIndex) {
    auto* devTraction1{this->traction1};
    auto* devTraction2{this->traction2};
    auto* devFaultStresses{this->faultStresses};
    auto* devTractionResults{this->tractionResults};

    sycl::nd_range rng{{this->currLayerSize * misc::NumPaddedPoints}, {misc::NumPaddedPoints}};
    this->queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
        const auto ltsFace = item.get_group().get_group_id(0);
        const auto pointIndex = item.get_local_id(0);

        auto& faultStresses = devFaultStresses[ltsFace];
        auto& tractionResults = devTractionResults[ltsFace];

        // calculate traction
        tractionResults.traction1[timeIndex][pointIndex] =
            faultStresses.traction1[timeIndex][pointIndex];
        tractionResults.traction2[timeIndex][pointIndex] =
            faultStresses.traction2[timeIndex][pointIndex];
        devTraction1[ltsFace][pointIndex] = tractionResults.traction1[timeIndex][pointIndex];
        devTraction2[ltsFace][pointIndex] = tractionResults.traction2[timeIndex][pointIndex];
      });
    });
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  void saveDynamicStressOutput() {}

  void preHook(real (*stateVariableBuffer)[misc::NumPaddedPoints]) {}
  void postHook(real (*stateVariableBuffer)[misc::NumPaddedPoints]) {}

  protected:
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_
