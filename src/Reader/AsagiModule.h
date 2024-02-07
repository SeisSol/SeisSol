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

#ifndef ASAGI_MODULE_H
#define ASAGI_MODULE_H

#ifdef USE_ASAGI

#include "Parallel/MPI.h"

#include <string>

#include <asagi.h>

#include "utils/env.h"
#include "utils/logger.h"

#include "Modules/Module.h"
#include "Modules/Modules.h"

namespace seissol::asagi {

enum MPI_Mode { MPI_OFF, MPI_WINDOWS, MPI_COMM_THREAD, MPI_UNKNOWN };

class AsagiModule : public Module {
  private:
  /** The MPI mode used for ASAGI communication */
  MPI_Mode m_mpiMode;

  /** The real name set via the environment variable */
  std::string m_mpiModeName;

  /** The total number of threads (including the communication thread */
  int m_totalThreads;
  AsagiModule();

  public:
  AsagiModule(const AsagiModule&) = delete;
  void operator=(const AsagiModule&) = delete;

  void preMPI() override;

  /**
   * At the moment this function will only be registered when this warning/error
   * needs to be emitted.
   */
  void postMPIInit() override;

  /**
   * This hook is only registered if the comm thread is required
   */
  void preModel() override;

  /**
   * This hook is only registered if the comm thread is required
   */
  void postModel() override;

  static AsagiModule& getInstance();

  private:
  /**
   * First guess on the MPI mode
   *
   * The final decision will be made before starting the communication thread
   *
   * @warning This function is called before MPI initialization
   */
  static MPI_Mode getMPIMode();

  /**
   * @warning This function is called before MPI initialization
   */
  static int getTotalThreads();

  public:
  /**
   * @return The MPI mode for ASAGI
   */
  static MPI_Mode mpiMode();

  /**
   * @return The total number of threads available for ASAGI
   */
  static int totalThreads();

  private:
  static const char* ENV_MPI_MODE;
};

} // namespace seissol::asagi

#endif // USE_ASAGI

#endif // ASAGI_MODULE_H
