// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITYRECORDER_H_
#define SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITYRECORDER_H_

#include "Core/Common/Real.h"
#include "Core/Kernels/BatchRecorders/Recorders.h"
#include "Core/Memory/Descriptor/LTS.h"

namespace seissol::recording {

/// Records the batches of the plasticity of a layer of the configuration `Cfg`.
template <typename Cfg>
class PlasticityRecorder : public AbstractRecorder {
  public:
  using real = Real<Cfg>;

  protected:
  void setUpContext(LTS::Layer& layer) override { AbstractRecorder::setUpContext(layer); }

  public:
  void record(LTS::Layer& layer) override;
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITYRECORDER_H_
