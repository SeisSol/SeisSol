/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
 **/

#include "NRFReader.h"
#include <utils/logger.h>

#include <netcdf.h>

#include <cassert>

void check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    logError() << "line" << line << "of" << file << ":" << nc_strerror(stat) << std::endl;
  }
}

void seissol::sourceterm::readNRF(char const* filename, NRF& nrf)
{
  int ncid;
  int stat;

  /* dimension ids */
  int source_dim;
  int sroffset_dim;
  int sample1_dim;
  int sample2_dim;
  int sample3_dim;

  /* dimension lengths */
  size_t sroffset_len;
  size_t sample_len[3];

  /* variable ids */
  int centres_id;
  int subfaults_id;
  int sroffsets_id;
  int sliprates1_id;
  int sliprates2_id;
  int sliprates3_id;

  /* open nrf */
  stat = nc_open(filename, NC_NOWRITE, &ncid);
  check_err(stat,__LINE__,__FILE__);

  /* get dimensions */
  stat = nc_inq_dimid(ncid, "source", &source_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_inq_dimlen(ncid, source_dim, &nrf.source);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_inq_dimid(ncid, "sroffset", &sroffset_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_inq_dimlen(ncid, sroffset_dim, &sroffset_len);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_inq_dimid(ncid, "sample1", &sample1_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_inq_dimlen(ncid, sample1_dim, &sample_len[0]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_inq_dimid(ncid, "sample2", &sample2_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_inq_dimlen(ncid, sample2_dim, &sample_len[1]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_inq_dimid(ncid, "sample3", &sample3_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_inq_dimlen(ncid, sample3_dim, &sample_len[2]);
  check_err(stat,__LINE__,__FILE__);

  assert( nrf.source + 1 == sroffset_len );

  /* get varids */
  stat = nc_inq_varid(ncid, "centres", &centres_id);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_inq_varid(ncid, "subfaults", &subfaults_id);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_inq_varid(ncid, "sroffsets", &sroffsets_id);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_inq_varid(ncid, "sliprates1", &sliprates1_id);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_inq_varid(ncid, "sliprates2", &sliprates2_id);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_inq_varid(ncid, "sliprates3", &sliprates3_id);
  check_err(stat,__LINE__,__FILE__);

  /* allocate memory */
  assert(sizeof(glm::dvec3) == 3*sizeof(double));
  nrf.centres = new glm::dvec3[nrf.source];
  nrf.sroffsets = new Offsets[nrf.source + 1];
  nrf.subfaults = new Subfault[nrf.source];
  for (unsigned i = 0; i < 3; ++i) {
    nrf.sliprates[0] = new double[sample_len[0]];
    nrf.sliprates[1] = new double[sample_len[1]];
    nrf.sliprates[2] = new double[sample_len[2]];
  }

  /* get values */
  stat = nc_get_var(ncid, centres_id, nrf.centres);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_get_var(ncid, sroffsets_id, nrf.sroffsets);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_get_var(ncid, subfaults_id, nrf.subfaults);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_get_var_double(ncid, sliprates1_id, nrf.sliprates[0]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_get_var_double(ncid, sliprates2_id, nrf.sliprates[1]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_get_var_double(ncid, sliprates3_id, nrf.sliprates[2]);
  check_err(stat,__LINE__,__FILE__);

  /* close nrf */
  stat = nc_close(ncid);
  check_err(stat,__LINE__,__FILE__);
}
