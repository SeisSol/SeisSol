// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_INTEGRATIONSTAGE_H_
#define SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_INTEGRATIONSTAGE_H_

#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/CellStage.h"

#include <cstddef>

namespace seissol::kernels {

/**
 * The time integral of the unknowns of the cells of a cell cluster, for the outputs of the wave
 * field that write it: each correction adds the unknowns it leaves times the width of its step
 * (only first order in time).
 */
template <typename Cfg>
class IntegrationStage : public seissol::solver::CellStage<Cfg> {
  public:
  explicit IntegrationStage(LTS::Layer& layer) : layer_(layer) {}

  void beginCorrection(const seissol::solver::StepParams& params) override;
  std::size_t correctCell(std::size_t cell, LTS::Ref<Cfg>& data) override;
  void correctDevice(const seissol::solver::StepParams& params,
                     bool recurring,
                     parallel::runtime::StreamRuntime& runtime) override;

  private:
  LTS::Layer& layer_;
  double timestep_{0};
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_INTEGRATIONSTAGE_H_
