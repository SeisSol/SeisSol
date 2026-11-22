// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PlasticityStage.h"

#include "Config.h"
#include "Core/Common/Marker.h"
#include "Core/Kernels/Common.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Features/Plasticity/Solver/Plasticity.h"
#include "Features/Plasticity/Solver/Storage.h"

#include <cstddef>

#ifdef ACL_DEVICE
#include "Core/Common/BasicTypedefs.h"
#include "Core/Parallel/DeviceGraph.h"
#endif

namespace seissol::kernels {

template <typename Cfg>
PlasticityStage<Cfg>::PlasticityStage(RunContext& run,
                                      double tv,
                                      CompoundGlobalData<Cfg> globalData,
                                      LTS::Layer& layer)
    : run_(run), tv_(tv), globalData_(globalData), layer_(layer),
      yieldedHost_(1, seissol::memory::Memkind::Standard),
      yieldedDevice_(1,
                     isDeviceOn() ? seissol::memory::Memkind::DeviceGlobalMemory
                                  : seissol::memory::Memkind::Standard) {
  const auto* cellInformation = layer_.var<LTS::CellInformation>();
  for (std::size_t i = 0; i < layer_.size(); ++i) {
    if (cellInformation[i].plasticityEnabled) {
      ++plasticCells_;
    }
  }

  const auto [check, yield] = Plasticity<Cfg>::metrics();
  check_ = check;
  yield_ = yield;
  checkHandle_ = run_.flopCounter().addMetric("plasticity-check", "PL");
  yieldHandle_ = run_.flopCounter().addMetric("plasticity-yield", "PL");

  yieldedHost_[0] = 0;
  yieldedDevice_.copyFrom(yieldedHost_);
}

template <typename Cfg>
void PlasticityStage<Cfg>::beginCorrection(const seissol::solver::StepParams& params) {
  timestep_ = params.timeStepSize;
  oneMinusIntegratingFactor_ = Plasticity<Cfg>::computeRelaxTime(tv_, timestep_);
  plasticity_ = layer_.var<LTSPlasticity::Plasticity>(Cfg());
  pstrain_ = layer_.var<LTSPlasticity::PStrain>(Cfg());
}

template <typename Cfg>
std::size_t PlasticityStage<Cfg>::correctCell(std::size_t cell, LTS::Ref<Cfg>& data) {
  if (!data.template get<LTS::CellInformation>().plasticityEnabled) {
    return 0;
  }
  return Plasticity<Cfg>::computePlasticity(oneMinusIntegratingFactor_,
                                            timestep_,
                                            tv_,
                                            globalData_.onHost,
                                            &plasticity_[cell],
                                            data.template get<LTS::Dofs>(),
                                            pstrain_[cell]);
}

template <typename Cfg>
void PlasticityStage<Cfg>::endCorrection(std::size_t count) {
  yieldedHost_[0] += count;
}

template <typename Cfg>
void PlasticityStage<Cfg>::correctDevice(
    SEISSOL_GPU_PARAM const seissol::solver::StepParams& params,
    SEISSOL_GPU_PARAM bool recurring,
    SEISSOL_GPU_PARAM parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  const double timeStepWidth = params.timeStepSize;
  auto& table = layer_.getConditionalTable();
  const auto plasticityGraphKey =
      initializer::GraphKey(ComputeGraphType::Plasticity, timeStepWidth);
  auto* plasticity =
      layer_.var<LTSPlasticity::Plasticity>(Cfg(), seissol::initializer::AllocationPlace::Device);
  auto* isAdjustableVector =
      layer_.var<LTSPlasticity::FlagScratch>(seissol::initializer::AllocationPlace::Device);
  runtime.runGraph(
      plasticityGraphKey,
      layer_,
      [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
        Plasticity<Cfg>::computePlasticityBatched(timeStepWidth,
                                                  tv_,
                                                  globalData_.onDevice,
                                                  table,
                                                  plasticity,
                                                  yieldedDevice_.data(),
                                                  isAdjustableVector,
                                                  streamRuntime);
      },
      recurring);
#endif
}

template <typename Cfg>
void PlasticityStage<Cfg>::countCorrection() {
  // here instead of with the device work, which a replay of a recording skips; the yielding cells
  // get counted on the device, and added at the end of the phase
  run_.flopCounter().incrementMetric(checkHandle_, check_ * plasticCells_);
}

template <typename Cfg>
void PlasticityStage<Cfg>::finishPhase() {
  run_.flopCounter().incrementMetric(yieldHandle_, yield_ * yieldedHost_[0]);

  yieldedHost_.copyFrom(yieldedDevice_);
  run_.flopCounter().incrementMetric(yieldHandle_, yield_ * yieldedHost_[0]);

  yieldedHost_[0] = 0;
  yieldedDevice_.copyFrom(yieldedHost_);
}

#define SEISSOL_INSTANTIATE(Cfg) template class PlasticityStage<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels
