/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 * Read a mesh file in memory efficient way
 **/

#ifndef MESH_READER_H
#define MESH_READER_H

#include "MeshDefinition.h"
#include "MeshTools.h"

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <mpi.h>
#include <vector>
#include <array>
#include <unordered_map>
#include "Parallel/MPI.h"

namespace seissol::geometry {

enum class MeshFormat : int { Netcdf, PUML };

struct TransferRegion {
  int rank;
  std::size_t start;
  std::size_t size;
};

struct PassThrough {};

class MeshReader {
  protected:
  const int m_rank;

  std::vector<Element> m_elements;

  std::vector<Vertex> m_vertices;

  std::vector<std::reference_wrapper<const BoundaryElement>> m_boundaryElements;

  std::vector<TransferRegion> m_regions;

  std::unordered_map<int, std::vector<BoundaryElement>> m_MPINeighbors;

  /** Convert global element index to local */
  std::map<int, int> m_g2lElements;

  /** Convert global vertex index to local */
  std::map<int, int> m_g2lVertices;

  /** Number of MPI fault neighbors */
  std::map<int, std::vector<MPINeighborElement>> m_MPIFaultNeighbors;

  /** Fault information */
  std::vector<Fault> m_fault;

  /** Has a plus fault side */
  bool m_hasPlusFault;

  protected:
  MeshReader(int rank);

  public:
  virtual ~MeshReader();

  const std::vector<Element>& getElements() const;
  const std::vector<Vertex>& getVertices() const;
  const std::unordered_map<int, std::vector<BoundaryElement>>& getMPINeighbors() const;
  const std::vector<std::reference_wrapper<const BoundaryElement>>& getBoundaryElements() const;
  const std::map<int, std::vector<MPINeighborElement>>& getMPIFaultNeighbors() const;
  const std::vector<Fault>& getFault() const;
  bool hasFault() const;
  bool hasPlusFault() const;

  void displaceMesh(const std::array<double, 3>& displacement);

  // scalingMatrix is stored column-major, i.e.
  // scalingMatrix_ij = scalingMatrix[j][i]
  void scaleMesh(const std::array<std::array<double, 3>, 3>& scalingMatrix);

  /**
   * Reconstruct the fault information from the boundary conditions
   */
  void extractFaultInformation(const VrtxCoords refPoint, const int refPointMethod);

  void exchangeGhostlayerMetadata();

  template <typename SendT, typename OutT, typename InT, typename PreF, typename PostF>
  std::vector<OutT> exchangeGhostData(const std::vector<InT>& data,
                                      MPI_Comm comm,
                                      MPI_Datatype datatype,
                                      PreF&& preSelector = PassThrough(),
                                      PostF&& postSelector = PassThrough()) {
#ifdef USE_MPI
    const std::size_t size = m_boundaryElements.size();
    std::vector<SendT> sendElements(size);
    std::vector<SendT> recvElements(size);

    if constexpr (std::is_same_v<PreF, PassThrough>) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < size; ++i) {
        sendElements[i] = data[m_boundaryElements[i].get().localElement];
      }
    } else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < size; ++i) {
        sendElements[i] =
            std::forward<SendT>(std::invoke(preSelector,
                                            i,
                                            m_boundaryElements[i].get().localElement,
                                            data[m_boundaryElements[i].get().localElement]));
      }
    }

    std::vector<MPI_Request> requests(m_regions.size() * 2);
    std::size_t index = 0;
    const int tag = 15;
    for (const auto& region : m_regions) {
      MPI_Isend(&sendElements[region.start],
                region.size,
                datatype,
                region.rank,
                tag,
                comm,
                &requests[index]);
      MPI_Irecv(&recvElements[region.start],
                region.size,
                datatype,
                region.rank,
                tag,
                comm,
                &requests[index + 1]);
      index += 2;
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

    if constexpr (std::is_same_v<PostF, PassThrough>) {
      return recvElements;
    } else {
      std::vector<OutT> outElements(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t i = 0; i < size; ++i) {
        outElements[i] = std::forward<OutT>(std::invoke(postSelector, recvElements[i]));
      }
      return outElements;
    }
#else
    return {};
#endif
  }
};

} // namespace seissol::geometry

#endif // MESH_READER_H
