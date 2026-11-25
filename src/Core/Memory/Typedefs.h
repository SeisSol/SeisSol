// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_CORE_MEMORY_TYPEDEFS_H_
#define SEISSOL_SRC_CORE_MEMORY_TYPEDEFS_H_

#include "Alignment.h"
#include "CellLocalInformation.h"
#include "Config.h"
#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/Constants.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Common/Real.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Solver/Data.h"
#include "WavePropagation/Solver/SolverSelector.h"

#include <cstddef>

namespace seissol {

namespace kernels {
template <typename Cfg>
constexpr std::size_t NumSpaceQuadraturePoints =
    (Cfg::ConvergenceOrder + 1) * (Cfg::ConvergenceOrder + 1);
} // namespace kernels

/**
 * The generated constant matrices of the configuration `Cfg`, as one table of pointers into the
 * pool of its generator.
 *
 * There are two of these per configuration: one built on the image in this binary, one on a
 * copy of it in device memory. Which entries a table holds is decided by the code generator.
 **/
template <typename Cfg>
using GlobalData = seissol::Pool<Cfg>;

/// The global data of the configuration `Cfg` on the host and, in a GPU build, on the device.
template <typename Cfg>
struct CompoundGlobalData {
  GlobalData<Cfg>* onHost{nullptr};
  GlobalData<Cfg>* onDevice{nullptr};
};

// data for the cell local integration
template <typename Cfg>
struct alignas(Alignment) LocalIntegrationData {
  // star matrices
  Real<Cfg> starMatrices[3][seissol::tensor::star<Cfg>::size(0)]{};

  // flux solver for element local contribution
  Real<Cfg> nApNm1[4][seissol::tensor::AplusT<Cfg>::size()]{};

  // solver-specific data
  seissol::kernels::SolverOf<Cfg>::LocalData specific;
};

// data for the neighboring boundary integration
template <typename Cfg>
struct alignas(Alignment) NeighboringIntegrationData {
  // flux solver for the contribution of the neighboring elements
  Real<Cfg> nAmNm1[4][seissol::tensor::AminusT<Cfg>::size()]{};

  // solver-specific data
  seissol::kernels::SolverOf<Cfg>::NeighborData specific;
};

template <typename Cfg>
struct CellDRMapping {
  std::int8_t side{};
  std::int8_t faceRelation{};
  Real<Cfg>* godunov{nullptr};
  Real<Cfg>* fluxSolver{nullptr};
};

// The tensors start aligned, as the generated kernels may read them: dirichletOffset, for one, with
// the fused simulations first, which a kernel of libxsmm loads aligned.
template <typename Cfg>
struct BoundaryFaceInformation {
  // nodes is an array of 3d-points in global coordinates.
  alignas(Alignment)
      Real<Cfg> nodes[seissol::nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisDim<Cfg>] * 3]{};
  alignas(Alignment) Real<Cfg> dataT[seissol::tensor::T<Cfg>::size()]{};
  alignas(Alignment) Real<Cfg> dataTinv[seissol::tensor::Tinv<Cfg>::size()]{};
  alignas(Alignment) Real<Cfg> dirichletOffset[seissol::tensor::dirichletOffset<Cfg>::size()]{};
  alignas(Alignment) Real<Cfg> dirichletMap[seissol::tensor::dirichletMap<Cfg>::size()]{};
  Real<Cfg> fsgData[3]{};
};

template <typename Cfg>
struct CellBoundaryMapping {
  Real<Cfg>* nodes{nullptr};
  Real<Cfg>* dataT{nullptr};
  Real<Cfg>* dataTinv{nullptr};
  Real<Cfg>* dirichletOffset{nullptr};
  Real<Cfg>* dirichletMap{nullptr};
  Real<Cfg>* fsgData{nullptr};

  CellBoundaryMapping() = default;
  explicit CellBoundaryMapping(BoundaryFaceInformation<Cfg>& faceInfo)
      : nodes(faceInfo.nodes), dataT(faceInfo.dataT), dataTinv(faceInfo.dataTinv),
        dirichletOffset(faceInfo.dirichletOffset), dirichletMap(faceInfo.dirichletMap),
        fsgData(faceInfo.fsgData) {}
};

} // namespace seissol

#endif // SEISSOL_SRC_CORE_MEMORY_TYPEDEFS_H_
