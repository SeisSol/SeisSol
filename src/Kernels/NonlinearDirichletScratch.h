// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARDIRICHLETSCRATCH_H_
#define SEISSOL_SRC_KERNELS_NONLINEARDIRICHLETSCRATCH_H_

#include "GeneratedCode/tensor.h"

#include <cstddef>

namespace seissol::kernels {

/// What the device kernels of the nonlinear Dirichlet boundary keep per face, between the time
/// derivatives of its cell and its ghost state (LTS::NonlinearDirichletScratch): the state inside
/// at the time points of a step, a tensor like INodal (or INodalAligned) per time point; the state
/// of the cell at one of them on the way there, like I; and the ghost state in the face-aligned
/// basis, like INodalAligned. In reals, from the start of the block of the face.
template <typename Cfg>
struct NonlinearDirichletScratch {
  static constexpr std::size_t Inner = 0;
  static constexpr std::size_t State = Cfg::ConvergenceOrder * tensor::INodal<Cfg>::size();
  static constexpr std::size_t Aligned = State + tensor::I<Cfg>::size();
  static constexpr std::size_t PerFace = Aligned + tensor::INodalAligned<Cfg>::size();
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_NONLINEARDIRICHLETSCRATCH_H_
