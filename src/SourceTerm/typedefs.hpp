/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
 * Copyright (c) 2023, Intel corporation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Point source computation.
 **/

#ifndef SOURCETERM_TYPEDEFS_HPP_
#define SOURCETERM_TYPEDEFS_HPP_

#include <Initializer/MemoryAllocator.h>
#include <Kernels/precision.hpp>
#include <array>
#include <cstdlib>
#include <generated_code/tensor.h>
#include <memory>
#include <vector>

#ifdef ACL_DEVICE
#include "Device/UsmAllocator.h"
#endif

namespace seissol::kernels {
class PointSourceCluster;
} // namespace seissol::kernels

namespace seissol::sourceterm {
struct PointSourceClusterPair {
  std::unique_ptr<kernels::PointSourceCluster> host{nullptr};
  std::unique_ptr<kernels::PointSourceCluster> device{nullptr};
};

/** Models point sources of the form
 *    S(xi, eta, zeta, t) := (1 / |J|) * S(t) * M * delta(xi-xi_s, eta-eta_s, zeta-zeta_s),
 * where S(t) : t -> \mathbb R is the moment time history,
 * M \in \mathbb R^{9} contains entries of the moment tensor,
 * and delta is the 3-dimensional dirac distribution.
 *
 * (The scaling factor (1 / |J|) is due to the coordinate transformation (x,y,z) -> (xi,eta,zeta).)
 **/
struct PointSources {
  constexpr static unsigned TensorSize =
      (tensor::momentFSRM::Size > 9) ? tensor::momentFSRM::Size : 9;
  enum Mode { NRF, FSRM };
  enum Mode mode = NRF;

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
  seissol::memory::MemkindArray<std::array<real, 81>> stiffnessTensor;

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
};
} // namespace seissol::sourceterm

#endif
