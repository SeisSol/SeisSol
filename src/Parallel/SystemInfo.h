// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_PARALLEL_SYSTEMINFO_H_
#define SEISSOL_SRC_PARALLEL_SYSTEMINFO_H_

#include <string>
#include <vector>

namespace seissol {

/**
 * Provides infos about the system we run on.
 * Also prints them on the command line on setup.
 * TODO: maybe move most printing code from SeisSol.h here.
 */
class SystemInfo {
  public:
  /**
   * Collects data.
   */
  void init();

  /**
   * @return hostnames for all ranks in the communicator of the application
   */
  [[nodiscard]] const auto& getHostNames() const { return hostNames_; }

  [[nodiscard]] const auto& getPCIAddresses() const { return pcis_; }

  private:
  void loadHostInfo();
  void loadCPUInfo();
  void loadGPUInfo();

  std::vector<std::string> hostNames_;
  std::vector<std::string> pcis_;
};

} // namespace seissol
#endif // SEISSOL_SRC_PARALLEL_SYSTEMINFO_H_
