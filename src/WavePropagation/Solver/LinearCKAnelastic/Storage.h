// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_LINEARCKANELASTIC_STORAGE_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_LINEARCKANELASTIC_STORAGE_H_

#include "Alignment.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Real.h"
#include "Core/Common/TensorTraits.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "GeneratedCode/tensor.h"

/**
 * The variables of the cells that the anelastic solver keeps in the storage of the cells (see
 * LTS): the anelastic unknowns, which it keeps apart from the elastic ones, and its scratchpads on
 * the device.
 */
namespace seissol::LTS {

// empty if the configuration has no Qane
template <typename Cfg>
using DofsAneArray = Real<Cfg>[zeroGuard(kernels::size<tensor::Qane<Cfg>>())];

struct DofsAne : public initializer::VariantVariable<DofsAneArray> {};
struct IDofsAneScratch : public initializer::VariantScratchpad<Real> {};
struct DofsExtScratch : public initializer::VariantScratchpad<Real> {};

} // namespace seissol::LTS

namespace seissol::kernels::solver::linearckanelastic {

/// Adds the variables of the anelastic solver to the storage of the cells `storage`. The anelastic
/// unknowns get memory only if a configuration has some, but are there in any case, since the
/// outputs read them wherever they exist.
inline void addStorage(LTS::Storage& storage) {
  using namespace initializer;
  bool anelastic = false;
  forEachConfig([&](auto cfg) { anelastic |= kernels::size<tensor::Qane<decltype(cfg)>>() > 0; });
  if (anelastic) {
    storage.add<LTS::DofsAne>(
        LayerMask(Ghost), PagesizeHeap, LTS::allocationModeWP(LTS::AllocationPreset::Dofs));
  } else {
    storage.add<LTS::DofsAne>(LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior),
                              PagesizeHeap,
                              LTS::allocationModeWP(LTS::AllocationPreset::Dofs));
  }

  if constexpr (isDeviceOn()) {
    storage.add<LTS::IDofsAneScratch>(LayerMask(), Alignment, AllocationMode::DeviceOnly);
    storage.add<LTS::DofsExtScratch>(LayerMask(), Alignment, AllocationMode::DeviceOnly);
  }
}

/// Sets the sizes of the scratchpads of the anelastic solver in the layer `layer` of cells of the
/// configuration `Cfg`.
template <typename Cfg>
void setScratchpadSizes(LTS::Layer& layer) {
  using real = Real<Cfg>;
  layer.setEntrySize<LTS::IDofsAneScratch>(layer.size() * kernels::size<tensor::Iane<Cfg>>() *
                                           sizeof(real));
  layer.setEntrySize<LTS::DofsExtScratch>(layer.size() * kernels::size<tensor::Qext<Cfg>>() *
                                          sizeof(real));
}

} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_LINEARCKANELASTIC_STORAGE_H_
