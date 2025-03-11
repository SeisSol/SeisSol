// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_GLOBALDATA_H_
#define SEISSOL_SRC_INITIALIZER_GLOBALDATA_H_

#include "Initializer/Typedefs.h"
#include "MemoryAllocator.h"
#include <yateto.h>

#ifdef ACL_DEVICE
#include "device.h"
#endif // ACL_DEVICE

namespace seissol::initializer {
/*
 * \class MemoryProperties
 *
 * \brief An auxiliary data structure for a policy-based design
 *
 * Attributes are initialized with CPU memory properties by default.
 * See, an example of a policy-based design in GlobalData.cpp
 * */
struct MemoryProperties {
  size_t alignment{Alignment};
  size_t pagesizeHeap{PagesizeHeap};
  size_t pagesizeStack{PagesizeStack};
};

namespace matrixmanip {
struct OnHost {
  using CopyManagerT = typename yateto::DefaultCopyManager<real>;
  static MemoryProperties getProperties();
  static void negateStiffnessMatrix(GlobalData& globalData);
  static void initSpecificGlobalData(GlobalData& globalData,
                                     memory::ManagedAllocator& allocator,
                                     CopyManagerT& copyManager,
                                     size_t alignment,
                                     seissol::memory::Memkind memkind);
};

struct OnDevice {
  struct DeviceCopyPolicy {
    static real* copy(const real* first, const real* last, real*& mem);
  };
  using CopyManagerT = typename yateto::CopyManager<real, DeviceCopyPolicy>;
  static MemoryProperties getProperties();
  static void negateStiffnessMatrix(GlobalData& globalData);
  static void initSpecificGlobalData(GlobalData& globalData,
                                     memory::ManagedAllocator& allocator,
                                     CopyManagerT& copyManager,
                                     size_t alignment,
                                     seissol::memory::Memkind memkind);
};
} // namespace matrixmanip

// Generalized Global data initializers of SeisSol.
template <typename MatrixManipPolicyT>
struct GlobalDataInitializer {
  static void init(GlobalData& globalData,
                   memory::ManagedAllocator& memoryAllocator,
                   enum memory::Memkind memkind);
};

// Specific Global data initializers of SeisSol.
using GlobalDataInitializerOnHost = GlobalDataInitializer<matrixmanip::OnHost>;
using GlobalDataInitializerOnDevice = GlobalDataInitializer<matrixmanip::OnDevice>;
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_GLOBALDATA_H_
