// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

#ifndef SEISSOL_SRC_READER_ASAGIMODULE_H_
#define SEISSOL_SRC_READER_ASAGIMODULE_H_

#ifdef USE_ASAGI

#include "Parallel/MPI.h"

#include <string>

#include <asagi.h>

#include "utils/env.h"
#include "utils/logger.h"

#include "Modules/Module.h"
#include "Modules/Modules.h"

namespace seissol::asagi {

enum class AsagiMPIMode { Off, Windows, CommThread, Unknown };

class AsagiModule : public Module {
  private:
  /** The MPI mode used for ASAGI communication */
  AsagiMPIMode m_mpiMode;

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
  static AsagiMPIMode getMPIMode();

  /**
   * @warning This function is called before MPI initialization
   */
  static int getTotalThreads();

  public:
  /**
   * @return The MPI mode for ASAGI
   */
  static AsagiMPIMode mpiMode();

  /**
   * @return The total number of threads available for ASAGI
   */
  static int totalThreads();

  private:
  static inline const char* EnvMPIMode = "SEISSOL_ASAGI_MPI_MODE";
};

} // namespace seissol::asagi

#endif // USE_ASAGI

#endif // SEISSOL_SRC_READER_ASAGIMODULE_H_
