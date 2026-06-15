// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_READER_ASAGIMODULE_H_
#define SEISSOL_SRC_READER_ASAGIMODULE_H_

#ifdef USE_ASAGI

#include "Modules/Module.h"
#include "Parallel/Pin.h"

#include <memory>
#include <string>
#include <utils/env.h>

namespace seissol::asagi {

enum class AsagiMPIMode { Off, Windows, CommThread, Unknown };

class AsagiModule : public Module {
  private:
  utils::Env env_;

  /** The MPI mode used for ASAGI communication */
  AsagiMPIMode mpiMode_;

  /** The real name set via the environment variable */
  std::string mpiModeName_;

  /** The total number of threads (including the communication thread */
  int totalThreads_;

  parallel::Pinning* pinning_{};

  static std::shared_ptr<AsagiModule> instance;

  public:
  explicit AsagiModule(utils::Env& env);
  ~AsagiModule() override = default;
  AsagiModule(const AsagiModule&) = delete;
  void operator=(const AsagiModule&) = delete;
  AsagiModule(AsagiModule&&) = delete;
  void operator=(AsagiModule&&) = delete;

  utils::Env& getEnv();

  void preMPI() override;

  /**
   * At the moment this function will only be registered when this warning/error
   * needs to be emitted.
   */
  void postMPIInit() override;

  /**
   * This hook is only registered if the comm thread is required
   */
  void preMesh() override;

  /**
   * This hook is only registered if the comm thread is required
   */
  void postModel() override;

  static void initInstance(utils::Env& env);

  static AsagiModule& getInstance();

  private:
  /**
   * First guess on the MPI mode
   *
   * The final decision will be made before starting the communication thread
   *
   * @warning This function is called before MPI initialization
   */
  static AsagiMPIMode getMPIMode(utils::Env& env);

  /**
   * @warning This function is called before MPI initialization
   */
  static int getTotalThreads(utils::Env& env);

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
  static inline const std::string EnvMpiMode = "ASAGI_MPI_MODE";
};

} // namespace seissol::asagi

#endif // USE_ASAGI

#endif // SEISSOL_SRC_READER_ASAGIMODULE_H_
