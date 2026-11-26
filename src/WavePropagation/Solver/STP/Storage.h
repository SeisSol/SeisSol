// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_STP_STORAGE_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_STP_STORAGE_H_

#include "Alignment.h"
#include "Core/Common/Real.h"
#include "Core/Common/TensorTraits.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "GeneratedCode/tensor.h"

/**
 * The variables of the cells that the space-time predictor keeps in the storage of the cells (see
 * LTS): the inverses of the Z matrices of its cells, as a scratchpad on the device.
 */
namespace seissol::LTS {

struct ZinvExtra : public initializer::VariantScratchpad<Real> {};

} // namespace seissol::LTS

namespace seissol::kernels::solver::stp {

/// Adds the variables of the space-time predictor to the storage of the cells `storage`.
inline void addStorage(LTS::Storage& storage) {
  if constexpr (isDeviceOn()) {
    storage.add<LTS::ZinvExtra>(
        initializer::LayerMask(), Alignment, initializer::AllocationMode::HostDevicePinned);
  }
}

/// Sets the sizes of the scratchpads of the space-time predictor in the layer `layer` of cells of
/// the configuration `Cfg`.
template <typename Cfg>
void setScratchpadSizes(LTS::Layer& layer) {
  layer.setEntrySize<LTS::ZinvExtra>(layer.size() * kernels::familySize<tensor::Zinv<Cfg>>() *
                                     sizeof(Real<Cfg>));
}

} // namespace seissol::kernels::solver::stp

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_STP_STORAGE_H_
