// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_GLOBALDATA_H_
#define SEISSOL_SRC_MEMORY_GLOBALDATA_H_

#include "Initializer/Typedefs.h"
#include "MemoryAllocator.h"

namespace seissol::initializer {

namespace matrixmanip {
/**
 * Where the pool is read from.
 *
 * On the host there is nothing to do: the image is in the binary already, and
 * a table built on it needs neither an allocation nor a copy. On the device it
 * is one allocation and one copy for the whole image, whatever it contains.
 * */
struct OnHost {
  static GlobalData pool(memory::ManagedAllocator& allocator, memory::Memkind memkind);
};

struct OnDevice {
  static GlobalData pool(memory::ManagedAllocator& allocator, memory::Memkind memkind);
};
} // namespace matrixmanip

template <typename MatrixManipPolicyT>
struct GlobalDataInitializer {
  static void init(GlobalData& globalData,
                   memory::ManagedAllocator& memoryAllocator,
                   enum memory::Memkind memkind);
};

using GlobalDataInitializerOnHost = GlobalDataInitializer<matrixmanip::OnHost>;
using GlobalDataInitializerOnDevice = GlobalDataInitializer<matrixmanip::OnDevice>;
} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_GLOBALDATA_H_
