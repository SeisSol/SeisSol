/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#ifndef GLOBALDATA_H_
#define GLOBALDATA_H_

#include "MemoryAllocator.h"
#include "Typedefs.h"
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
    real* copy(const real* first, const real* last, real*& mem);
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

#endif
