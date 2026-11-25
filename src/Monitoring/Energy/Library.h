// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_ENERGY_LIBRARY_H_
#define SEISSOL_SRC_MONITORING_ENERGY_LIBRARY_H_

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace seissol::monitoring::energy {

/**
 * @brief A shared library opened at runtime, which stays loaded for the rest of the run.
 *
 * The management libraries of the GPU vendors come with their drivers, not with SeisSol, and do not
 * depend on the runtime SeisSol is built for: SYCL on an NVIDIA GPU reads it through NVML as well.
 * A machine without them runs without their sensors.
 */
class Library {
  public:
  //! @brief Opens the first of @p names that the dynamic loader finds.
  static std::optional<Library> open(const std::vector<std::string>& names);

  //! @brief The function named @p name, which has the type @p Function; null if there is none.
  template <typename Function>
  [[nodiscard]] Function* function(const char* name) const {
    // dlsym gives functions as data pointers, which POSIX allows to convert back
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    return reinterpret_cast<Function*>(symbol(name));
  }

  //! @brief The name it was opened under.
  [[nodiscard]] const std::string& name() const { return name_; }

  private:
  Library(void* handle, std::string name) : handle_(handle), name_(std::move(name)) {}

  [[nodiscard]] void* symbol(const char* name) const;

  void* handle_;
  std::string name_;
};

/**
 * @brief The names to look for the library @p name of ROCm under: as it is, then in the library
 * directories of ROCM_PATH and of /opt/rocm.
 */
std::vector<std::string> rocmLibraryNames(const std::string& name);

} // namespace seissol::monitoring::energy

#endif // SEISSOL_SRC_MONITORING_ENERGY_LIBRARY_H_
