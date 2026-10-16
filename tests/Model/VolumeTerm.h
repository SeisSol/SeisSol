// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_VOLUMETERM_H_
#define SEISSOL_TESTS_MODEL_VOLUMETERM_H_

#include "Alignment.h"
#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Interface.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <yateto.h>

namespace seissol::unit_test {

namespace volumeterm {
/// Whether a kernel takes the volume term the predictor of a cell forms for its corrector
/// (LocalTmp::volumeFused): the predictor writes it, the corrector reads it.
template <typename KernelT, typename = void>
struct TakesVolumeFused : std::false_type {};

template <typename KernelT>
struct TakesVolumeFused<KernelT, std::void_t<decltype(std::declval<KernelT>().volumeFused)>>
    : std::true_type {};
} // namespace volumeterm

/// Hands a predictor or a corrector the buffer the volume term goes through, where the build
/// forms the volume term in the predictor.
template <typename KernelT>
void bindVolumeFused(KernelT& krnl, kernels::LocalTmp& tmp) {
  if constexpr (volumeterm::TakesVolumeFused<KernelT>::value) {
    krnl.volumeFused = tmp.volumeFused;
  }
}

/// Adds the volume term of the corrector for the time integral `timeIntegrated` to `update`, with
/// the operator `local` carries, the way a cell forms it.
///
/// Where the build forms the volume term in the predictor, the predictor runs first, with the
/// first Taylor weight one and the others zero: its time integral is then the field it starts
/// from, so the volume term it hands on is the one of `timeIntegrated`.
template <typename VolumeT, typename LocalIntegrationT>
void addVolumeTerm(real* timeIntegrated, real* update, const LocalIntegrationT& local) {
  constexpr bool FromPredictor = volumeterm::TakesVolumeFused<VolumeT>::value;

  VolumeT volume{};
  volume.bindGlobals(seissol::Pool::host());
  volume.Q = update;
  kernels::LocalTmp tmp(0.0);
  if constexpr (FromPredictor) {
    // named through the volume kernel, so that a build whose predictor hands nothing on never
    // instantiates a body that names what it would hand on
    using DerivativeT = std::conditional_t<FromPredictor, kernel::derivative, void>;
    alignas(Alignment) std::array<real, yateto::computeFamilySize<tensor::dQ>()> derivatives{};
    alignas(Alignment) std::array<real, tensor::I::size()> integrated{};

    DerivativeT derivative{};
    derivative.bindGlobals(seissol::Pool::host());
    kernels::bindStarOperands(derivative, local);
    kernels::bindSourceOperands(derivative, local);
    derivative.dQ(0) = timeIntegrated;
    for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      derivative.dQ(i) = derivatives.data() + yateto::computeFamilySize<tensor::dQ>(1, i);
    }
    derivative.I = integrated.data();
    for (std::size_t i = 0; i < ConvergenceOrder; ++i) {
      derivative.power(i) = i == 0 ? 1 : 0;
    }
    bindVolumeFused(derivative, tmp);
    derivative.execute();
  } else {
    volume.I = timeIntegrated;
    kernels::bindStarOperands(volume, local);
    kernels::bindSourceOperands(volume, local);
  }
  bindVolumeFused(volume, tmp);
  volume.execute();
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_VOLUMETERM_H_
