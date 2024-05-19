/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
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
 * Velocity field reader Fortran interface
 */

#ifndef ASAGIREADER_H
#define ASAGIREADER_H

#ifdef USE_ASAGI

#include "Parallel/MPI.h"
#include "easi/util/AsagiReader.h"
#include <asagi.h>

#include "utils/env.h"
#include "utils/logger.h"

#include "AsagiModule.h"
#include "Monitoring/instrumentation.hpp"

namespace seissol::asagi {
enum class NumaCacheMode { Off, On, Cache };

class AsagiReader : public easi::AsagiReader {
  private:
  /** Prefix for environment variables */
  const std::string envPrefix;

  /** Number of threads used by ASAGI */
  unsigned int asagiThreads;

#ifdef USE_MPI
  /** MPI communicator used by ASAGI */
  MPI_Comm comm;
#endif

  public:
  AsagiReader(const char* envPrefix
#ifdef USE_MPI
              ,
              MPI_Comm comm = seissol::MPI::mpi.comm()
#endif
  );

  virtual ::asagi::Grid* open(const char* file, const char* varname);
  virtual unsigned numberOfThreads() const;

  private:
  static NumaCacheMode getNumaMode();
};

} // namespace seissol::asagi

#endif

#endif // ASAGIREADER_H
