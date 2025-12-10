// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_GLOBALDATA_H_
#define SEISSOL_SRC_MEMORY_GLOBALDATA_H_

#include "Config.h"
#include "Initializer/Typedefs.h"
#include "MemoryAllocator.h"

#include <Common/Executor.h>
#include <Common/Templating.h>
#include <variant>
#include <yateto.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif // ACL_DEVICE

namespace seissol {

struct GlobalData {
  private:
  memory::ManagedAllocator allocator;
  enum memory::Memkind memkindHost = memory::Memkind::HighBandwidth; // TODO: really?
  enum memory::Memkind memkindDevice = memory::Memkind::DeviceGlobalMemory;

  using GlobalDataArray = ChangeVariadicT<
      std::tuple,
      TransformVariadicT<std::optional, TransformVariadicT<GlobalDataCfg, ConfigVariant>>>;

  GlobalDataArray onHost;
  GlobalDataArray onDevice;

  public:
  GlobalData() = default;

  void init(std::size_t configId);

  template <typename Cfg, Executor Exec = Executor::Host>
  const GlobalDataCfg<Cfg>& get() const {
    if constexpr (Exec == Executor::Host) {
      return std::get<std::optional<GlobalDataCfg<Cfg>>>(onHost).value();
    } else if constexpr (Exec == Executor::Device) {
      return std::get<std::optional<GlobalDataCfg<Cfg>>>(onDevice).value();
    }
    throw;
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_GLOBALDATA_H_
