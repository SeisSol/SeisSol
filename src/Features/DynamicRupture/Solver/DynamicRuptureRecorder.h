// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURERECORDER_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURERECORDER_H_

#include "Core/Common/Real.h"
#include "Core/Kernels/BatchRecorders/Recorders.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <unordered_map>

namespace seissol::recording {

/// Records the batches of a layer of dynamic rupture faces of the configuration `Cfg`.
template <typename Cfg>
class DynamicRuptureRecorder : public AbstractRecorder {
  public:
  using real = Real<Cfg>;

  void record(DynamicRupture::Layer& layer) override;

  protected:
  void setUpContext(DynamicRupture::Layer& layer) override {
    AbstractRecorder::setUpContext(layer);
  }

  private:
  void recordSpaceInterpolation();
  std::unordered_map<real*, real*> idofsAddressRegistry_;
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURERECORDER_H_
