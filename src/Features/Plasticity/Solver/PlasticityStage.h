// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITYSTAGE_H_
#define SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITYSTAGE_H_

#include "Core/Common/Real.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/FlopCounter.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/CellStage.h"
#include "Features/Plasticity/Physics/Plasticity.h"
#include "Features/Plasticity/Solver/Storage.h"

#include <cstddef>

namespace seissol::kernels {

/**
 * The plasticity of the cells of a cell cluster, which yield after each correction where the
 * stress exceeds the yield criterion; with the relaxation time `tv`.
 */
template <typename Cfg>
class PlasticityStage : public seissol::solver::CellStage<Cfg> {
  public:
  using real = Real<Cfg>;

  PlasticityStage(monitoring::FlopCounter& flopCounter,
                  double tv,
                  CompoundGlobalData<Cfg> globalData,
                  LTS::Layer& layer);

  void beginCorrection(const seissol::solver::StepParams& params) override;
  std::size_t correctCell(std::size_t cell, LTS::Ref<Cfg>& data) override;
  void endCorrection(std::size_t count) override;
  void correctDevice(const seissol::solver::StepParams& params,
                     bool recurring,
                     parallel::runtime::StreamRuntime& runtime) override;
  void countCorrection() override;
  void finishPhase() override;

  private:
  monitoring::FlopCounter& flopCounter_;
  double tv_;
  CompoundGlobalData<Cfg> globalData_;
  LTS::Layer& layer_;

  // what the cells share in the current correction
  double timestep_{0};
  real oneMinusIntegratingFactor_{0};
  const model::PlasticityData<Cfg>* plasticity_{nullptr};
  LTSPlasticity::PStrainArray<Cfg>* pstrain_{nullptr};

  // the cells with plasticity, the FLOPs of the check of a cell and of the yielding of a cell
  std::size_t plasticCells_{0};
  PerformanceEstimate check_{};
  PerformanceEstimate yield_{};
  std::size_t checkHandle_{};
  std::size_t yieldHandle_{};

  // the cells that yielded since the last phase, on the host and on the device
  memory::MemkindArray<std::size_t> yieldedHost_;
  memory::MemkindArray<std::size_t> yieldedDevice_;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITYSTAGE_H_
