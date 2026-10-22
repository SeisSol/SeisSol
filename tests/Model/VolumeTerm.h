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

/// Whether a kernel reads the source term as the one matrix a cell carries in the data of its
/// solver (model::initializeSpecificLocalData), rather than at the sample points of its material.
template <typename KernelT, typename = void>
struct ReadsCellSource : std::false_type {};

template <typename KernelT>
struct ReadsCellSource<KernelT, std::void_t<decltype(std::declval<KernelT>().ET)>>
    : std::true_type {};

/// Whether the volume term of any time integral can be formed the way a cell forms it: by the
/// volume kernel, or by the derivative chain run on that time integral as the field (see
/// addVolumeTerm). A predictor that solves for its space-time field hands on the volume term of
/// that field alone.
template <typename VolumeT>
constexpr bool FormsVolumeTermOf =
    !TakesVolumeFused<VolumeT>::value || TakesVolumeFused<kernel::derivative>::value;
} // namespace volumeterm

/// Hands a predictor or a corrector the buffer the volume term goes through, where the build
/// forms the volume term in the predictor.
template <typename KernelT>
void bindVolumeFused(KernelT& krnl, kernels::LocalTmp& tmp) {
  if constexpr (volumeterm::TakesVolumeFused<KernelT>::value) {
    krnl.volumeFused = tmp.volumeFused;
  }
}

/// Hands a kernel the source term of the cell in whichever shape the cell carries it, as the host
/// does: the one matrix of the cell, or the scalars at the sample points of its material.
template <typename KernelT, typename LocalIntegrationT>
void bindSourceTerm(KernelT& krnl, const LocalIntegrationT& local) {
  if constexpr (volumeterm::ReadsCellSource<KernelT>::value) {
    krnl.ET = local.specific.sourceMatrix;
  }
  kernels::bindSourceOperands(krnl, local);
}

/// Adds the volume term of the corrector for the time integral `timeIntegrated` to `update`, with
/// the operator `local` carries, the way a cell forms it.
///
/// Where the build forms the volume term in the predictor, the predictor runs first, with the
/// first Taylor weight one and the others zero: its time integral is then the field it starts
/// from, so the volume term it hands on is the one of `timeIntegrated`.
template <typename VolumeT, typename LocalIntegrationT>
void addVolumeTerm(real* timeIntegrated, real* update, const LocalIntegrationT& local) {
  static_assert(volumeterm::FormsVolumeTermOf<VolumeT>,
                "the predictor of this build forms the volume term of its own field alone");
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
    bindSourceTerm(derivative, local);
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
    bindSourceTerm(volume, local);
  }
  bindVolumeFused(volume, tmp);
  volume.execute();
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_VOLUMETERM_H_
