// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_NEIGHINTEGRATIONRECORDER_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_NEIGHINTEGRATIONRECORDER_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Real.h"
#include "Core/Kernels/BatchRecorders/Recorders.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "GeneratedCode/configboundary.h"
#include "WavePropagation/Physics/Common.h"

#include <array>
#include <cstddef>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::recording {

/// Records the batches of the neighbor integration of a layer of the configuration `Cfg`.
template <typename Cfg>
class NeighIntegrationRecorder : public AbstractRecorder {
  public:
  using real = Real<Cfg>;

  void record(LTS::Layer& layer) override;

  protected:
  void setUpContext(LTS::Layer& layer) override {
    integratedDofsAddressCounter_ = 0;
    configBoundaryBatches_.assign(builtConfigCount(), {});
    fromCanonical_ = {};
    converted_ = {};
    fromCoupledCanonical_ = {};
    coupledConverted_ = {};
    coupling_ = {};
    canonicalRegistry_.clear();
    convertedRegistry_.clear();
    configBoundaryScratchCounter_ = 0;
    AbstractRecorder::setUpContext(layer);
  }

  private:
  void recordDofsTimeEvaluation();
  void recordNeighborFluxIntegrals();

  /// Records the face `face` of the cell `cell`, whose neighbor computes in another configuration
  /// (see kernels::ConfigBoundary): the time integral of the neighbor in its configuration if it
  /// provides derivatives, its conversion into the canonical form, and from there into the
  /// configuration of the layer for the side of the neighbor; each once per neighbor (and side).
  void recordConfigBoundaryFace(std::size_t cell, std::size_t face);
  /// Sets the batches recordConfigBoundaryFace collected.
  void recordConfigBoundaryBatches();
  /// The next `bytes` of the scratch of the faces between configurations.
  void* allocateConfigBoundaryScratch(std::size_t bytes);

  std::unordered_map<real*, real*> idofsAddressRegistry_;
  size_t integratedDofsAddressCounter_{0};

  // the canonical forms, in the same precision for all configurations
  using Canonical = generated::ConfigBoundaryKernels<Cfg>::CanonicalReal;

  // per configuration of the neighbors
  struct ConfigBoundaryBatches {
    // the derivatives of the neighbors and their time integrals, with the GTS and the LTS relation
    std::array<std::vector<void*>, 2> derivatives;
    std::array<std::vector<void*>, 2> integrals;
    // the time integrals of the neighbors to convert, and their canonical forms
    std::vector<void*> toCanonical;
    std::vector<Canonical*> canonical;
  };
  std::vector<ConfigBoundaryBatches> configBoundaryBatches_;
  // per side of the neighbors: the canonical forms, and their conversions
  std::array<std::vector<Canonical*>, Cell::NumFaces> fromCanonical_;
  std::array<std::vector<real*>, Cell::NumFaces> converted_;
  // the same per family of the neighbors of other families (see model::familyIndex), with the
  // coupling matrices of the faces
  template <typename T>
  using PerFamilyAndSide =
      std::array<std::array<std::vector<T>, Cell::NumFaces>, model::FamilyCount>;
  PerFamilyAndSide<Canonical*> fromCoupledCanonical_;
  PerFamilyAndSide<real*> coupledConverted_;
  PerFamilyAndSide<real*> coupling_;
  // the canonical form of each neighbor, and its conversion per side
  std::unordered_map<const void*, Canonical*> canonicalRegistry_;
  std::map<std::pair<const void*, std::size_t>, real*> convertedRegistry_;
  std::size_t configBoundaryScratchCounter_{0};
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_NEIGHINTEGRATIONRECORDER_H_
