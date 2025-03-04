// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_SOURCETERM_TYPEDEFS_H_
#define SEISSOL_SRC_SOURCETERM_TYPEDEFS_H_

#include "Common/Constants.h"
#include "Kernels/Precision.h"
#include "Memory/MemoryAllocator.h"
#include "generated_code/tensor.h"
#include <array>
#include <cstdlib>

namespace seissol::sourceterm {
/** Models point sources of the form
 *    S(xi, eta, zeta, t) := (1 / |J|) * S(t) * M * delta(xi-xi_s, eta-eta_s, zeta-zeta_s),
 * where S(t) : t -> \mathbb R is the moment time history,
 * M \in \mathbb R^{9} contains entries of the moment tensor,
 * and delta is the 3-dimensional dirac distribution.
 *
 * (The scaling factor (1 / |J|) is due to the coordinate transformation (x,y,z) -> (xi,eta,zeta).)
 **/

enum class PointSourceMode { Nrf, Fsrm };

struct PointSources {
  constexpr static unsigned TensorSize =
      (tensor::momentFSRM::Size > 9) ? tensor::momentFSRM::Size : 9;
  PointSourceMode mode = PointSourceMode::Nrf;

  /** mInvJInvPhisAtSources[][k] := M_{kl}^-1 * |J|^-1 * phi_l(xi_s, eta_s, zeta_s), where phi_l is
   * the l-th basis function and xi_s, eta_s, and zeta_s are the space position
   *  of the point source in the reference tetrahedron. */
  seissol::memory::MemkindArray<
      seissol::memory::AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>>
      mInvJInvPhisAtSources;

  /** NRF: Basis vectors of the fault.
   * 0-2: Tan1X-Z   = first fault tangent (main slip direction in most cases)
   * 3-5: Tan2X-Z   = second fault tangent
   * 6-8: NormalX-Z = fault normal
   *
   * FSRM: Moment tensor */
  seissol::memory::MemkindArray<seissol::memory::AlignedArray<real, TensorSize>> tensor;

  /// Area
  seissol::memory::MemkindArray<real> A;

  /// elasticity tensor
  seissol::memory::MemkindArray<seissol::memory::AlignedArray<real, 81>> stiffnessTensor;

  /// onset time
  seissol::memory::MemkindArray<double> onsetTime;

  /// sampling interval
  seissol::memory::MemkindArray<double> samplingInterval;

  /// offset into slip rate vector
  std::array<seissol::memory::MemkindArray<std::size_t>, 3u> sampleOffsets;

  /** NRF: slip rate in
   * 0: Tan1 direction
   * 1: Tan2 direction
   * 2: Normal direction
   *
   * FSRM: 0: slip rate (all directions) */
  std::array<seissol::memory::MemkindArray<real>, 3u> sample;

  /** Number of point sources in this struct. */
  unsigned numberOfSources = 0;

  PointSources(seissol::memory::Memkind memkind)
      : mInvJInvPhisAtSources(memkind), tensor(memkind), A(memkind), stiffnessTensor(memkind),
        onsetTime(memkind), samplingInterval(memkind),
        sampleOffsets{seissol::memory::MemkindArray<std::size_t>(memkind),
                      seissol::memory::MemkindArray<std::size_t>(memkind),
                      seissol::memory::MemkindArray<std::size_t>(memkind)},
        sample{seissol::memory::MemkindArray<real>(memkind),
               seissol::memory::MemkindArray<real>(memkind),
               seissol::memory::MemkindArray<real>(memkind)} {}
  PointSources(const PointSources& source, seissol::memory::Memkind memkind)
      : mInvJInvPhisAtSources(source.mInvJInvPhisAtSources, memkind),
        tensor(source.tensor, memkind), A(source.A, memkind),
        stiffnessTensor(source.stiffnessTensor, memkind), onsetTime(source.onsetTime, memkind),
        samplingInterval(source.samplingInterval, memkind),
        sampleOffsets{seissol::memory::MemkindArray<std::size_t>(source.sampleOffsets[0], memkind),
                      seissol::memory::MemkindArray<std::size_t>(source.sampleOffsets[1], memkind),
                      seissol::memory::MemkindArray<std::size_t>(source.sampleOffsets[2], memkind)},
        sample{seissol::memory::MemkindArray<real>(source.sample[0], memkind),
               seissol::memory::MemkindArray<real>(source.sample[1], memkind),
               seissol::memory::MemkindArray<real>(source.sample[2], memkind)} {}
  ~PointSources() { numberOfSources = 0; }
};

struct CellToPointSourcesMapping {
  //! Pointer to DOFs
  real (*dofs)[tensor::Q::size()];
  //! First point source that has an effect on the cell
  unsigned pointSourcesOffset;
  /** The point sources buffer is ordered by cells, hence the point sources
   * that affect the cell with copyInteriorOffset reside in
   * {pointSourcesOffset, ..., pointSourcesOffset + numberOfPointSources - 1}
   * in the point sources buffer.
   **/
  unsigned numberOfPointSources;

  CellToPointSourcesMapping() : dofs(nullptr), pointSourcesOffset(0), numberOfPointSources(0) {}
};

struct ClusterMapping {
  seissol::memory::MemkindArray<unsigned> sources;
  seissol::memory::MemkindArray<CellToPointSourcesMapping> cellToSources;

  ClusterMapping(seissol::memory::Memkind memkind) : sources(memkind), cellToSources(memkind) {}
  ClusterMapping(const ClusterMapping& mapping, seissol::memory::Memkind memkind)
      : sources(mapping.sources, memkind), cellToSources(mapping.cellToSources, memkind) {}
};
} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_TYPEDEFS_H_
