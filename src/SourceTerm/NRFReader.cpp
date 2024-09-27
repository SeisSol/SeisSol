/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
 * Copyright (c) 2023, Intel Corporation
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
 **/

#include "NRFReader.h"
#include <SourceTerm/NRF.h>
#include <cstddef>
#include <ostream>
#include <utils/logger.h>

#include <netcdf.h>

#include <cassert>

namespace {

static void check_err(const int stat, const int line, const char* file) {
  if (stat != NC_NOERR) {
    logError() << "line" << line << "of" << file << ":" << nc_strerror(stat) << std::endl;
  }
}

} // namespace

void seissol::sourceterm::readNRF(const char* filename, NRF& nrf) {
  int ncid;
  int stat;

  /* dimension ids */
  int sourceDim;
  int sroffsetDim;
  int sample1Dim;
  int sample2Dim;
  int sample3Dim;

  /* dimension lengths */
  size_t sourceLen;
  size_t sroffsetLen;
  size_t sampleLen[3];

  /* variable ids */
  int centresId;
  int subfaultsId;
  int sroffsetsId;
  int sliprates1Id;
  int sliprates2Id;
  int sliprates3Id;

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
