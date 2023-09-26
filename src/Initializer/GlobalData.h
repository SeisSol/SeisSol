/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
#include "typedefs.hpp"
#include <Common/configs.hpp>
#include <Common/templating.hpp>
#include <Common/varianthelper.hpp>
#include <optional>
#include <yateto.h>

#include "Common/executor.hpp"

#ifdef ACL_DEVICE
#include "device.h"
#endif // ACL_DEVICE

namespace seissol {
  namespace initializers {
    namespace matrixmanip {
      template<typename Config>
      struct OnHost {
        using RealT = typename Config::RealT;
        using CopyManagerT = typename yateto::DefaultCopyManager<RealT>;
        static MemoryProperties getProperties();
        static void negateStiffnessMatrix(GlobalData<Config>& globalData);
        static void initSpecificGlobalData(GlobalData<Config>& globalData,
                                           memory::ManagedAllocator& allocator,
                                           CopyManagerT& copyManager,
                                           size_t alignment,
                                           seissol::memory::Memkind memkind);
      };

      template<typename Config>
      struct OnDevice {
        using RealT = typename Config::RealT;
        struct DeviceCopyPolicy {
          RealT* copy(RealT const* first, RealT const* last, RealT*& mem);
        };
        using CopyManagerT = typename yateto::CopyManager<RealT, DeviceCopyPolicy>;
        static MemoryProperties getProperties();
        static void negateStiffnessMatrix(GlobalData<Config>& globalData);
        static void initSpecificGlobalData(GlobalData<Config>& globalData,
                                           memory::ManagedAllocator& allocator,
                                           CopyManagerT& copyManager,
                                           size_t alignment,
                                           seissol::memory::Memkind memkind);
      };
    }  // namespace matrixmanip

    // Generalized Global data initializers of SeisSol.
    template<typename Config, typename MatrixManipPolicyT>
    struct GlobalDataInitializer {
      using RealT = typename Config::RealT;
      static void init(GlobalData<Config> &globalData,
                       memory::ManagedAllocator &memoryAllocator,
                       enum memory::Memkind memkind);
    };

    // Specific Global data initializers of SeisSol.
    template<typename Config>
    using GlobalDataInitializerOnHost = GlobalDataInitializer<Config, matrixmanip::OnHost<Config>>;

    template<typename Config>
    using GlobalDataInitializerOnDevice = GlobalDataInitializer<Config, matrixmanip::OnDevice<Config>>;

    class GlobalDataStorage {
    public:
      using GlobalDataArray = ChangeVariadicT<std::tuple, TransformVariadicT<std::optional, TransformVariadicT<GlobalData, SupportedConfigs>>>;

      GlobalDataStorage(memory::ManagedAllocator& allocator) : allocator(allocator) {}

      template<typename Executor, typename Config>
      const GlobalData<Config>& getData() {
        constexpr auto typeId = variantTypeId<Config, SupportedConfigs>();
        if constexpr (Executor == Executor::Host) {
          auto& element = std::get<typeId>(onHost);
          if (!std::get<typeId>(onHost).has_value()) {
            GlobalData<Config> globalData;
            GlobalDataInitializerOnHost<Config>::init(globalData, allocator, memkindHost);
            element = std::move(globalData);
          }
          return element;
        }
        else if constexpr (Executor == Executor::Device) {
          auto& element = std::get<typeId>(onHost);
          if (!std::get<typeId>(onDevice).has_value()) {
            GlobalData<Config> globalData;
            GlobalDataInitializerOnDevice<Config>::init(globalData, allocator, memkindDevice);
            element = std::move(globalData);
          }
          return element;
        }
      }
    private:
      memory::ManagedAllocator& allocator;
      enum memory::Memkind memkindHost = memory::Memkind::HighBandwidth; // TODO: really?
      enum memory::Memkind memkindDevice = memory::Memkind::DeviceGlobalMemory;

      GlobalDataArray onHost;
      GlobalDataArray onDevice;
    };
  }  // namespace initializers
} // namespace seissol

#endif
