// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_KERNELS_GLOBALDATASTORE_H_
#define SEISSOL_SRC_CORE_KERNELS_GLOBALDATASTORE_H_

#include "Core/Common/ConfigDispatch.h"
#include "Core/Kernels/Common.h"
#include "Core/Memory/Typedefs.h"

namespace seissol {

/// The global data of every configuration built into the executable: on the host and, in a GPU
/// build, on the device as well.
class GlobalDataStore {
  public:
  /// The global data of the configuration `Cfg` on both host and device.
  template <typename Cfg>
  [[nodiscard]] CompoundGlobalData<Cfg> get() {
    CompoundGlobalData<Cfg> global{};
    global.onHost = &onHost_.get<Cfg>();
    if constexpr (isDeviceOn()) {
      global.onDevice = &onDevice_.get<Cfg>();
    }
    return global;
  }

  template <typename Cfg>
  [[nodiscard]] GlobalData<Cfg>& onHost() {
    return onHost_.get<Cfg>();
  }

  template <typename Cfg>
  [[nodiscard]] GlobalData<Cfg>& onDevice() {
    return onDevice_.get<Cfg>();
  }

  private:
  PerConfig<GlobalData> onHost_;
  PerConfig<GlobalData> onDevice_;
};

} // namespace seissol

#endif // SEISSOL_SRC_CORE_KERNELS_GLOBALDATASTORE_H_
