/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 */

#ifndef MPI_BASIC_H
#define MPI_BASIC_H

namespace seissol {

/**
 * Basic MPI abstraction
 */
class MPIBasic {
  protected:
  /** This rank */
  int m_rank{0};

  /** Rank in the shared memory sub-communicator */
  int m_sharedMemMpiRank{0};

  /** Number of processors */
  int m_size{1};

  /** Number of ranks in the shared memory sub-communicator */
  int m_sharedMemMpiSize{1};

  MPIBasic() = default;

  public:
  virtual ~MPIBasic() = default;

  /**
   * @return The rank of this process
   */
  [[nodiscard]] int rank() const { return m_rank; }

  /**
   * @return The rank within the shared memory sub-communicator
   */
  [[nodiscard]] int sharedMemMpiRank() const { return m_sharedMemMpiRank; }

  /**
   * @return The total number of processes
   */
  [[nodiscard]] int size() const { return m_size; }

  /**
   * @return The number of ranks within the shared memory sub-communicator
   */
  [[nodiscard]] int sharedMemMpiSize() const { return m_sharedMemMpiSize; }

  [[nodiscard]] bool isSingleProcess() const { return size() == 1; }

  [[nodiscard]] bool isSingleNode() const { return size() == sharedMemMpiSize(); }
};
} // namespace seissol

#endif // MPI_BASIC_H
