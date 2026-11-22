// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "IntegrationStage.h"

#include "Config.h"
#include "Core/Common/Marker.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "GeneratedCode/tensor.h"

#include <cstddef>

#ifdef ACL_DEVICE
#include "Core/Common/Real.h"
#include "Core/Kernels/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Core/Kernels/BatchRecorders/DataTypes/EncodedConstants.h"

#include <Device/device.h>
#endif

namespace seissol::kernels {

template <typename Cfg>
void IntegrationStage<Cfg>::beginCorrection(const seissol::solver::StepParams& params) {
  timestep_ = params.timeStepSize;
}

template <typename Cfg>
std::size_t IntegrationStage<Cfg>::correctCell(std::size_t /*cell*/, LTS::Ref<Cfg>& data) {
  auto* __restrict integral = data.template get<LTS::Integrals>();
  const auto* __restrict dofs = data.template get<LTS::Dofs>();
  const auto timestep = timestep_;

// only first-order time integration for the output here
#pragma omp simd
  for (std::size_t dof = 0; dof < tensor::Q<Cfg>::size(); ++dof) {
    integral[dof] += timestep * dofs[dof];
  }
  return 0;
}

template <typename Cfg>
void IntegrationStage<Cfg>::correctDevice(
    const seissol::solver::StepParams& /*params*/,
    bool /*recurring*/,
    SEISSOL_GPU_PARAM parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;
  using real = Real<Cfg>;
  auto& table = layer_.getConditionalTable();
  const ConditionalKey key = ConditionalKey(*KernelNames::Time);
  if (table.contains(key)) {
    auto entry = table.at(key);
    device::DeviceInstance::instance().algorithms().accumulateBatchedData(
        const_cast<const real**>(entry.get<real*>(inner_keys::Wp::Id::Dofs)->getDeviceDataPtr()),
        entry.get<real*>(inner_keys::Wp::Id::Integrals)->getDeviceDataPtr(),
        tensor::Q<Cfg>::Size,
        entry.get<real*>(inner_keys::Wp::Id::Dofs)->getSize(),
        runtime.stream());
  }
#endif
}

#define SEISSOL_INSTANTIATE(Cfg) template class IntegrationStage<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels
