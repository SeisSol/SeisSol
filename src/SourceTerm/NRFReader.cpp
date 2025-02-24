// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "NRFReader.h"
#include <SourceTerm/NRF.h>
#include <cstddef>
#include <ostream>
#include <utils/logger.h>

#include <netcdf.h>

#include <cassert>

namespace {

void check_err(const int stat, const int line, const char* file) {
  if (stat != NC_NOERR) {
    logError() << "line" << line << "of" << file << ":" << nc_strerror(stat) << std::endl;
  }
}

} // namespace

void seissol::sourceterm::readNRF(const char* filename, NRF& nrf) {
  int ncid = 0;
  int stat = 0;

  /* dimension ids */
  int sourceDim = 0;
  int sroffsetDim = 0;
  int sample1Dim = 0;
  int sample2Dim = 0;
  int sample3Dim = 0;

  /* dimension lengths */
  size_t sourceLen = 0;
  size_t sroffsetLen = 0;
  size_t sampleLen[3];

  /* variable ids */
  int centresId = 0;
  int subfaultsId = 0;
  int sroffsetsId = 0;
  int sliprates1Id = 0;
  int sliprates2Id = 0;
  int sliprates3Id = 0;

  /* open nrf */
  stat = nc_open(filename, NC_NOWRITE, &ncid);
  check_err(stat, __LINE__, __FILE__);

  /* get dimensions */
  stat = nc_inq_dimid(ncid, "source", &sourceDim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_inq_dimlen(ncid, sourceDim, &sourceLen);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_dimid(ncid, "sroffset", &sroffsetDim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_inq_dimlen(ncid, sroffsetDim, &sroffsetLen);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_dimid(ncid, "sample1", &sample1Dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_inq_dimlen(ncid, sample1Dim, &sampleLen[0]);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_dimid(ncid, "sample2", &sample2Dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_inq_dimlen(ncid, sample2Dim, &sampleLen[1]);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_dimid(ncid, "sample3", &sample3Dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_inq_dimlen(ncid, sample3Dim, &sampleLen[2]);
  check_err(stat, __LINE__, __FILE__);

  assert(sourceLen + 1 == sroffsetLen);

  /* get varids */
  stat = nc_inq_varid(ncid, "centres", &centresId);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_varid(ncid, "subfaults", &subfaultsId);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_varid(ncid, "sroffsets", &sroffsetsId);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_varid(ncid, "sliprates1", &sliprates1Id);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_varid(ncid, "sliprates2", &sliprates2Id);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_varid(ncid, "sliprates3", &sliprates3Id);
  check_err(stat, __LINE__, __FILE__);

  /* allocate memory */
  static_assert(sizeof(Eigen::Vector3d) == 3 * sizeof(double),
                "sizeof(Eigen::Vector3d) does not equal 3*sizeof(double).");
  nrf.centres.resize(sourceLen);
  nrf.sroffsets.resize(sourceLen + 1);
  nrf.subfaults.resize(sourceLen);
  for (unsigned i = 0; i < 3; ++i) {
    nrf.sliprates[i].resize(sampleLen[i]);
  }

  /* get values */
  stat = nc_get_var(ncid, centresId, nrf.centres.data());
  check_err(stat, __LINE__, __FILE__);

  stat = nc_get_var(ncid, sroffsetsId, nrf.sroffsets.data());
  check_err(stat, __LINE__, __FILE__);

  stat = nc_get_var(ncid, subfaultsId, nrf.subfaults.data());
  check_err(stat, __LINE__, __FILE__);

  stat = nc_get_var_double(ncid, sliprates1Id, nrf.sliprates[0].data());
  check_err(stat, __LINE__, __FILE__);

  stat = nc_get_var_double(ncid, sliprates2Id, nrf.sliprates[1].data());
  check_err(stat, __LINE__, __FILE__);

  stat = nc_get_var_double(ncid, sliprates3Id, nrf.sliprates[2].data());
  check_err(stat, __LINE__, __FILE__);

  /* close nrf */
  stat = nc_close(ncid);
  check_err(stat, __LINE__, __FILE__);
}
