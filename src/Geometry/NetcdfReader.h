/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian rettenberger AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2013-2023, SeisSol Group
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
 * Read a mesh in the Netcdf format
 **/

#ifndef NETCDF_READER_H
#define NETCDF_READER_H

#ifdef USE_NETCDF

#include "MeshReader.h"

namespace seissol::geometry {

class NetcdfReader : public seissol::geometry::MeshReader {
  public:
  NetcdfReader(int rank, int nProcs, const char* meshFile);

  private:
  /**
   * Adds the information about an MPI neighbor
   *
   * @localID The local id of the neighborhood
   * @bndRank The neighbor rank
   * @elemSize Number of boundary elements
   * @bndElemLocalIds List of boundary ids
   */
  void addMPINeighbor(int localID, int bndRank, int elemSize, const int* bndElemLocalIds);

  /**
   * Finds all locals elements for each vertex
   */
  void findElementsPerVertex();

  /**
   * Switch to collective access for a netCDf variable
   */
  static void collectiveAccess(int ncFile, int ncVar);

  static void checkNcError(int error);
};

} // namespace seissol::geometry

#endif // USE_NETCDF

#endif // NETCDF_READER_H
