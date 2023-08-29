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

#include "GlobalData.h"
#include <generated_code/init.h>
#include <yateto.h>
#include <type_traits>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace seissol::initializers {
  namespace matrixmanip {
    MemoryProperties OnHost::getProperties() {
      // returns MemoryProperties initialized with default values i.e., CPU memory properties
      return MemoryProperties();
    }

    void OnHost::negateStiffnessMatrix(GlobalData<Config> &globalData) {
      for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
        RealT *matrix = const_cast<RealT *>(globalData.stiffnessMatricesTransposed(transposedStiffness));
        for (unsigned i = 0; i < seissol::Yateto<Config>::Init::kDivMT::size(transposedStiffness); ++i) {
          matrix[i] *= -1.0;
        }
      }
    }

    void OnHost::initSpecificGlobalData(GlobalData<Config>& globalData,
                                        memory::ManagedAllocator& allocator,
                                        CopyManagerT& copyManager,
                                        size_t alignment,
                                        seissol::memory::Memkind memkind) {
      // thread-local LTS integration buffers
      int l_numberOfThreads = 1;
#ifdef _OPENMP
        l_numberOfThreads = omp_get_max_threads();
#endif
        RealT *integrationBufferLTS
          = (RealT *) allocator.allocateMemory(l_numberOfThreads * (4 * Yateto<Config>::Tensor::I::size()) * sizeof(RealT),
                                                          alignment,
                                                          memkind);

      // initialize w.r.t. NUMA
      #ifdef _OPENMP
      #pragma omp parallel
      {
        size_t l_threadOffset = omp_get_thread_num() * (4 * Yateto<Config>::Tensor::I::size());
      #else
        size_t l_threadOffset = 0;
      #endif
        for (unsigned int l_dof = 0; l_dof < (4 * Yateto<Config>::Tensor::I::size()); l_dof++) {
          integrationBufferLTS[l_dof + l_threadOffset] = (RealT) 0.0;
        }
      #ifdef _OPENMP
      }
      #endif

      globalData.integrationBufferLTS = integrationBufferLTS;
    }

    MemoryProperties OnDevice::getProperties() {
      MemoryProperties prop{};
#ifdef ACL_DEVICE
      device::DeviceInstance& device = device::DeviceInstance::getInstance();
      prop.alignment = device.api->getGlobMemAlignment();
      prop.pagesizeHeap = prop.alignment;
      prop.pagesizeStack = prop.alignment;
#endif
      return prop;
    }

    void OnDevice::negateStiffnessMatrix(GlobalData<Config> &globalData) {
#ifdef ACL_DEVICE
      device::DeviceInstance& device = device::DeviceInstance::getInstance();
      for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
        const RealT scaleFactor = -1.0;
        device.algorithms.scaleArray(const_cast<RealT*>(globalData.stiffnessMatricesTransposed(transposedStiffness)),
                                     scaleFactor,
                                     Yateto<Config>::Init::kDivMT::size(transposedStiffness),
                                     device.api->getDefaultStream());
      }
#endif // ACL_DEVICE
    }
    void OnDevice::initSpecificGlobalData(GlobalData<Config>& globalData,
                                          memory::ManagedAllocator& allocator,
                                          CopyManagerT& copyManager,
                                          size_t alignment,
                                          seissol::memory::Memkind memkind) {
#ifdef ACL_DEVICE
      const size_t size = yateto::alignedUpper(Yateto<Config>::Tensor::replicateInitialLoadingM::size(),
                                               yateto::alignedReals<RealT>(alignment));
      RealT* plasticityStressReplication =
          static_cast<RealT*>(allocator.allocateMemory(size * sizeof(RealT),
                                                      alignment,
                                                      memkind));

      copyManager.template copyTensorToMemAndSetPtr<Yateto<Config>::Init::replicateInitialLoadingM>(plasticityStressReplication,
                                                                                    globalData.replicateStresses,
                                                                                    alignment);
#endif // ACL_DEVICE
    }

    RealT* OnDevice::DeviceCopyPolicy::copy(RealT const* first, RealT const* last, RealT*& mem) {
#ifdef ACL_DEVICE
      device::DeviceInstance& device = device::DeviceInstance::getInstance();
      const unsigned bytes = (last - first) * sizeof(RealT);
      device.api->copyTo(mem, first, bytes);
      mem += (last - first);
      return mem;
#else // ACL_DEVICE
      return nullptr;
#endif
    }

  } // namespace matrixmanip


template<typename MatrixManipPolicyT>
void GlobalDataInitializer<MatrixManipPolicyT>::init(GlobalData<Config>& globalData,
                                                     memory::ManagedAllocator& memoryAllocator,
                                                     enum seissol::memory::Memkind memkind) {
  MemoryProperties prop = MatrixManipPolicyT::getProperties();

  // We ensure that global matrices always start at an aligned memory address,
  // such that mixed cases with aligned and non-aligned global matrices do also work.
  unsigned globalMatrixMemSize = 0;
  globalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::kDivM>(yateto::alignedReals<RealT>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::kDivMT>(yateto::alignedReals<RealT>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::rDivM>(yateto::alignedReals<RealT>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::rT>(yateto::alignedReals<RealT>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::fMrT>(yateto::alignedReals<RealT>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::fP>(yateto::alignedReals<RealT>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::nodal::V3mTo2nFace>(yateto::alignedReals<RealT>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::project2nFaceTo3m>(yateto::alignedReals<RealT>(prop.alignment));

  globalMatrixMemSize += yateto::alignedUpper(Yateto<Config>::Tensor::evalAtQP::size(),  yateto::alignedReals<RealT>(prop.alignment));
  globalMatrixMemSize += yateto::alignedUpper(Yateto<Config>::Tensor::projectQP::size(), yateto::alignedReals<RealT>(prop.alignment));

  RealT* globalMatrixMem = static_cast<RealT*>(memoryAllocator.allocateMemory(globalMatrixMemSize * sizeof(RealT),
                                                                            prop.pagesizeHeap,
                                                                            memkind));

  RealT* globalMatrixMemPtr = globalMatrixMem;
  typename MatrixManipPolicyT::CopyManagerT copyManager;
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::kDivMT>(globalMatrixMemPtr, globalData.stiffnessMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::kDivM>(globalMatrixMemPtr, globalData.stiffnessMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::rDivM>(globalMatrixMemPtr, globalData.changeOfBasisMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::rT>(globalMatrixMemPtr, globalData.neighbourChangeOfBasisMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::fMrT>(globalMatrixMemPtr, globalData.localChangeOfBasisMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::fP>(globalMatrixMemPtr, globalData.neighbourFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::nodal::V3mTo2nFace>(globalMatrixMemPtr, globalData.V3mTo2nFace, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::project2nFaceTo3m>(globalMatrixMemPtr, globalData.project2nFaceTo3m, prop.alignment);

  copyManager.template copyTensorToMemAndSetPtr<Yateto<Config>::Init::evalAtQP>(globalMatrixMemPtr, globalData.evalAtQPMatrix, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<Yateto<Config>::Init::projectQP>(globalMatrixMemPtr, globalData.projectQPMatrix, prop.alignment);

  assert(globalMatrixMemPtr == globalMatrixMem + globalMatrixMemSize);

  // @TODO Integrate this step into the code generator
  MatrixManipPolicyT::negateStiffnessMatrix(globalData);

  // Dynamic Rupture global matrices
  unsigned drGlobalMatrixMemSize = 0;
  drGlobalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::V3mTo2nTWDivM>(yateto::alignedReals<RealT>(prop.alignment));
  drGlobalMatrixMemSize += yateto::computeFamilySize<typename Yateto<Config>::Init::V3mTo2n>(yateto::alignedReals<RealT>(prop.alignment));

  RealT* drGlobalMatrixMem = static_cast<RealT*>(memoryAllocator.allocateMemory(drGlobalMatrixMemSize  * sizeof(RealT),
                                                                              prop.pagesizeHeap,
                                                                              memkind));

  RealT* drGlobalMatrixMemPtr = drGlobalMatrixMem;
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::V3mTo2nTWDivM>(drGlobalMatrixMemPtr, globalData.nodalFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<typename Yateto<Config>::Init::V3mTo2n>(drGlobalMatrixMemPtr, globalData.faceToNodalMatrices, prop.alignment);

  assert(drGlobalMatrixMemPtr == drGlobalMatrixMem + drGlobalMatrixMemSize);

  // Plasticity global matrices
  unsigned plasticityGlobalMatrixMemSize = 0;
  plasticityGlobalMatrixMemSize += yateto::alignedUpper(Yateto<Config>::Tensor::v::size(),    yateto::alignedReals<RealT>(prop.alignment));
  plasticityGlobalMatrixMemSize += yateto::alignedUpper(Yateto<Config>::Tensor::vInv::size(), yateto::alignedReals<RealT>(prop.alignment));

  RealT* plasticityGlobalMatrixMem
      = static_cast<RealT*>(memoryAllocator.allocateMemory(plasticityGlobalMatrixMemSize * sizeof(RealT),
                                                          prop.pagesizeHeap,
                                                          memkind));

  RealT* plasticityGlobalMatrixMemPtr = plasticityGlobalMatrixMem;
  copyManager.template copyTensorToMemAndSetPtr<Yateto<Config>::Init::v>(plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrix, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<Yateto<Config>::Init::vInv>(plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrixInverse, prop.alignment);

  assert(plasticityGlobalMatrixMemPtr == plasticityGlobalMatrixMem + plasticityGlobalMatrixMemSize);

  MatrixManipPolicyT::initSpecificGlobalData(globalData,
                                             memoryAllocator,
                                             copyManager,
                                             prop.pagesizeStack,
                                             memkind);
}

template void GlobalDataInitializer<matrixmanip::OnHost>::init(GlobalData& globalData,
                                                               memory::ManagedAllocator& memoryAllocator,
                                                               enum memory::Memkind memkind);

template void GlobalDataInitializer<matrixmanip::OnDevice>::init(GlobalData& globalData,
                                                                 memory::ManagedAllocator& memoryAllocator,
                                                                 enum memory::Memkind memkind);

} // namespace seissol::initializers
