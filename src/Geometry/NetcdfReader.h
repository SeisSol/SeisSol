// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian rettenberger AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

#ifndef SEISSOL_SRC_GEOMETRY_NETCDFREADER_H_
#define SEISSOL_SRC_GEOMETRY_NETCDFREADER_H_

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

  private:
  /**
   * Switch to collective access for a netCDf variable
   */
  static void collectiveAccess(int ncFile, int ncVar);

  static void checkNcError(int error);
};

} // namespace seissol::geometry

#endif // USE_NETCDF

#endif // SEISSOL_SRC_GEOMETRY_NETCDFREADER_H_
