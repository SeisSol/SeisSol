// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_LOCALINTEGRATIONRECORDER_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_LOCALINTEGRATIONRECORDER_H_

#include "Core/Common/Real.h"
#include "Core/Kernels/BatchRecorders/Recorders.h"
#include "Core/Memory/Descriptor/LTS.h"

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace seissol::recording {

/// Records the batches of the local integration of a layer of the configuration `Cfg`.
template <typename Cfg>
class LocalIntegrationRecorder : public AbstractRecorder {
  public:
  using real = Real<Cfg>;

  explicit LocalIntegrationRecorder(double g) : g_(g) {}

  void record(LTS::Layer& layer) override;

  protected:
  double g_{9.81};

  void setUpContext(LTS::Layer& layer) override {
    integratedDofsAddressCounter_ = 0;
    derivativesAddressCounter_ = 0;
    AbstractRecorder::setUpContext(layer);
  }

  private:
  void recordTimeAndVolumeIntegrals();
  void recordFreeSurfaceGravityBc();
  void recordDirichletBc();
  void recordAnalyticalBc(LTS::Layer& layer);
  void recordLocalFluxIntegral();
  void recordDisplacements();

  std::unordered_map<size_t, real*> idofsAddressRegistry_;
  std::vector<real*> dQPtrs_;

  size_t integratedDofsAddressCounter_{0};
  size_t derivativesAddressCounter_{0};
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_LOCALINTEGRATIONRECORDER_H_
