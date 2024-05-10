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

#include "Kernels/precision.hpp"
#include "generated_code/tensor.h"
#include <array>
#include <vector>
#include "Common/constants.hpp"
#include "Kernels/precision.hpp"
#include "generated_code/tensor.h"
#include <cstdint>
#include <cstdlib>
#include <vector>

#ifdef ACL_DEVICE
#include "Device/UsmAllocator.h"
#endif

namespace seissol::sourceterm {
#ifdef ACL_DEVICE
using AllocatorT = device::UsmAllocator<real>;
#else
using AllocatorT = std::allocator<real>;
#endif
template <typename T>
using VectorT =
    std::vector<T, typename std::allocator_traits<AllocatorT>::template rebind_alloc<T>>;

template <typename T, std::size_t N>
class AlignedArray {
  public:
  inline T* data() { return data_; }
  inline const T* data() const { return data_; }
  constexpr T& operator[](std::size_t pos) { return data_[pos]; }
  constexpr const T& operator[](std::size_t pos) const { return data_[pos]; }
  constexpr std::size_t size() const noexcept { return N; }

  private:
  alignas(Alignment) T data_[N];
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
  VectorT<AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>> mInvJInvPhisAtSources;

  /** NRF: Basis vectors of the fault.
   * 0-2: Tan1X-Z   = first fault tangent (main slip direction in most cases)
   * 3-5: Tan2X-Z   = second fault tangent
   * 6-8: NormalX-Z = fault normal
   *
   * FSRM: Moment tensor */
  VectorT<AlignedArray<real, TensorSize>> tensor;

  /// Area
  VectorT<real> A;

  /// elasticity tensor
  VectorT<std::array<real, 81>> stiffnessTensor;

  /// onset time
  VectorT<double> onsetTime;

  /// sampling interval
  VectorT<double> samplingInterval;

  /// offset into slip rate vector
  std::array<VectorT<std::size_t>, 3u> sampleOffsets;

  /** NRF: slip rate in
   * 0: Tan1 direction
   * 1: Tan2 direction
   * 2: Normal direction
   *
   * FSRM: 0: slip rate (all directions) */
  std::array<VectorT<real>, 3u> sample;

  /** Number of point sources in this struct. */
  unsigned numberOfSources = 0;

  PointSources(const AllocatorT& alloc)
      : mInvJInvPhisAtSources(decltype(mInvJInvPhisAtSources)::allocator_type(alloc)),
        tensor(decltype(tensor)::allocator_type(alloc)), A(alloc),
        stiffnessTensor(decltype(stiffnessTensor)::allocator_type(alloc)),
        onsetTime(decltype(onsetTime)::allocator_type(alloc)),
        samplingInterval(decltype(samplingInterval)::allocator_type(alloc)),
        sampleOffsets{VectorT<std::size_t>(VectorT<std::size_t>::allocator_type(alloc)),
                      VectorT<std::size_t>(VectorT<std::size_t>::allocator_type(alloc)),
                      VectorT<std::size_t>(VectorT<std::size_t>::allocator_type(alloc))},
        sample{VectorT<real>(alloc), VectorT<real>(alloc), VectorT<real>(alloc)} {}
  ~PointSources() { numberOfSources = 0; }
};

struct CellToPointSourcesMapping {
  //! Pointer to DOFs
  real (*dofs)[tensor::Q::size()];
  //! First point source that has an effect on the cell
  unsigned pointSourcesOffset;
  /** The point sources buffer is ordered by cells, hence the point sources
   * that affect the cell with copyInteriorOffset reside in
   * [pointSourcesOffset, pointSourcesOffset + numberOfPointSources)
   * in the point sources buffer.
   **/
  unsigned numberOfPointSources;

  CellToPointSourcesMapping() : dofs(nullptr), pointSourcesOffset(0), numberOfPointSources(0) {}
};

struct ClusterMapping {
  VectorT<unsigned> sources;
  VectorT<CellToPointSourcesMapping> cellToSources;

  ClusterMapping(const AllocatorT& alloc)
      : sources(alloc), cellToSources(decltype(cellToSources)::allocator_type(alloc)) {}
};
} // namespace seissol::sourceterm

#endif
