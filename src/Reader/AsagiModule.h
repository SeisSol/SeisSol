// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_READER_ASAGIMODULE_H_
#define SEISSOL_SRC_READER_ASAGIMODULE_H_

#ifdef USE_ASAGI

#include <string>

#include <asagi.h>

#include "Modules/Module.h"

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
  ~AsagiModule() override = default;
  AsagiModule(const AsagiModule&) = delete;
  void operator=(const AsagiModule&) = delete;
  AsagiModule(AsagiModule&&) = delete;
  void operator=(AsagiModule&&) = delete;

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
  static inline const std::string EnvMpiMode = "SEISSOL_ASAGI_MPI_MODE";
};

} // namespace seissol::asagi

#endif // USE_ASAGI

#endif // SEISSOL_SRC_READER_ASAGIMODULE_H_
