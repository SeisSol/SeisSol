// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "GlobalData.h"

#include "Alignment.h"
#include "Config.h"
#include "GeneratedCode/pool.h"
#include "Initializer/Typedefs.h"
#include "Memory/MemoryAllocator.h"

#include <cstddef>

namespace seissol::initializer {

namespace matrixmanip {
template <typename Cfg>
GlobalData<Cfg> OnHost::pool(memory::ManagedAllocator& /*allocator*/, memory::Memkind /*memkind*/) {
  return seissol::Pool<Cfg>::host();
}

template <typename Cfg>
GlobalData<Cfg> OnDevice::pool(memory::ManagedAllocator& allocator, memory::Memkind memkind) {
  const std::size_t bytes = seissol::poolBytes<Cfg>();
  void* image = allocator.allocateMemory(bytes, PagesizeHeap, memkind);
  seissol::memory::memcopy(
      image, seissol::poolData<Cfg>(), bytes, memkind, memory::Memkind::Standard);
  return seissol::Pool<Cfg>::create(image);
}
} // namespace matrixmanip

template <typename MatrixManipPolicyT>
template <typename Cfg>
void GlobalDataInitializer<MatrixManipPolicyT>::init(GlobalData<Cfg>& globalData,
                                                     memory::ManagedAllocator& memoryAllocator,
                                                     enum seissol::memory::Memkind memkind) {
  globalData = MatrixManipPolicyT::template pool<Cfg>(memoryAllocator, memkind);
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template void GlobalDataInitializer<matrixmanip::OnHost>::init<Cfg>(                             \
      GlobalData<Cfg> & globalData,                                                                \
      memory::ManagedAllocator & memoryAllocator,                                                  \
      enum memory::Memkind memkind);                                                               \
  template void GlobalDataInitializer<matrixmanip::OnDevice>::init<Cfg>(                           \
      GlobalData<Cfg> & globalData,                                                                \
      memory::ManagedAllocator & memoryAllocator,                                                  \
      enum memory::Memkind memkind);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::initializer
