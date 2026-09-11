// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "GlobalData.h"

#include "Alignment.h"
#include "GeneratedCode/pool.h"
#include "Initializer/Typedefs.h"
#include "Memory/MemoryAllocator.h"

#include <cstddef>

namespace seissol::initializer {
namespace matrixmanip {

GlobalData OnHost::pool(memory::ManagedAllocator& /*allocator*/, memory::Memkind /*memkind*/) {
  return seissol::Pool::host();
}

GlobalData OnDevice::pool(memory::ManagedAllocator& allocator, memory::Memkind memkind) {
  const std::size_t bytes = seissol::poolBytes();
  void* image = allocator.allocateMemory(bytes, PagesizeHeap, memkind);
  seissol::memory::memcopy(image, seissol::poolData(), bytes, memkind, memory::Memkind::Standard);
  return seissol::Pool::create(image);
}

} // namespace matrixmanip

template <typename MatrixManipPolicyT>
void GlobalDataInitializer<MatrixManipPolicyT>::init(GlobalData& globalData,
                                                     memory::ManagedAllocator& memoryAllocator,
                                                     enum seissol::memory::Memkind memkind) {
  globalData = MatrixManipPolicyT::pool(memoryAllocator, memkind);
}

template void
    GlobalDataInitializer<matrixmanip::OnHost>::init(GlobalData& globalData,
                                                     memory::ManagedAllocator& memoryAllocator,
                                                     enum memory::Memkind memkind);

template void
    GlobalDataInitializer<matrixmanip::OnDevice>::init(GlobalData& globalData,
                                                       memory::ManagedAllocator& memoryAllocator,
                                                       enum memory::Memkind memkind);

} // namespace seissol::initializer
