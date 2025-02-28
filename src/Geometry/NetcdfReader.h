// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_NETCDFREADER_H_
#define SEISSOL_SRC_GEOMETRY_NETCDFREADER_H_

#ifdef USE_NETCDF

#include "MeshReader.h"

namespace seissol::geometry {

class NetcdfReader : public seissol::geometry::MeshReader {
  public:
  NetcdfReader(int rank, int nProcs, const char* meshFile);

  bool inlineTimestepCompute() const override;
  bool inlineClusterCompute() const override;

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

#endif // SEISSOL_SRC_GEOMETRY_NETCDFREADER_H_
