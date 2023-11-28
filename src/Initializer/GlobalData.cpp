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
#include "SeisSol.h"

#ifdef _OPENMP
#  include <omp.h>
#include <Kernels/Filter.h>
#include <memory>
#endif

namespace init = seissol::init;

namespace seissol::initializers {
  namespace matrixmanip {
    MemoryProperties OnHost::getProperties() {
      // returns MemoryProperties initialized with default values i.e., CPU memory properties
      return MemoryProperties();
    }

    void OnHost::negateStiffnessMatrix(GlobalData &globalData) {
      for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
        real *matrix = const_cast<real *>(globalData.stiffnessMatricesTransposed(transposedStiffness));
        for (unsigned i = 0; i < init::kDivMT::size(transposedStiffness); ++i) {
          matrix[i] *= -1.0;
        }
      }
    }

    void OnHost::initSpecificGlobalData(GlobalData& globalData,
                                        memory::ManagedAllocator& allocator,
                                        CopyManagerT& copyManager,
                                        size_t alignment,
                                        seissol::memory::Memkind memkind) {
      // thread-local LTS integration buffers
      int l_numberOfThreads = 1;
#ifdef _OPENMP
        l_numberOfThreads = omp_get_max_threads();
#endif
        real *integrationBufferLTS
          = (real *) allocator.allocateMemory(l_numberOfThreads * (4 * tensor::I::size()) * sizeof(real),
                                                          alignment,
                                                          memkind);

      // initialize w.r.t. NUMA
      #ifdef _OPENMP
      #pragma omp parallel
      {
        size_t l_threadOffset = omp_get_thread_num() * (4 * tensor::I::size());
      #else
        size_t l_threadOffset = 0;
      #endif
        for (unsigned int l_dof = 0; l_dof < (4 * tensor::I::size()); l_dof++) {
          integrationBufferLTS[l_dof + l_threadOffset] = (real) 0.0;
        }
      #ifdef _OPENMP
      }
      #endif

      globalData.integrationBufferLTS = integrationBufferLTS;
    }


     void OnHost::applyFilterMatrixToGlobalData(GlobalData& globalData, const seissol::kernels::Filter* filter) {
      //TODO(Lukas) Implement
        // TODO Add all matrices? What about DR? Point source? Viscoelastic? Poroelasticity source? Initial field?
        // localFluxNodal = lambda i: self.Q['kp'] <= self.Q['kp'] + self.db.project2nFaceTo3m[i]['kn'] * self.INodal['no'] * self.AminusT['op']

       //localFlux = lambda i: self.Q['kp'] <= self.Q['kp'] + self.db.rDivM[i][self.t('km')] * self.db.fMrT[i][self.t('ml')] * self.I['lq'] * self.AplusT['qp']
       // neighborFlux = lambda h, j, i: self.Q['kp'] <= self.Q['kp'] + self.db.rDivM[i][self.t('km')] * self.db.fP[h][self.t('mn')] * self.db.rT[j][self.t('nl')] * self.I['lq'] * self.AminusT['qp']

       //globalData.project2nFaceTo3m
       //m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;

      auto changeOfBasisMatricesView = std::array<yateto::DenseTensorView<2, real>, 4>({
        init::rDivM::view<0>::create(const_cast<real*>(globalData.changeOfBasisMatrices(0))),
        init::rDivM::view<1>::create(const_cast<real*>(globalData.changeOfBasisMatrices(1))),
        init::rDivM::view<2>::create(const_cast<real*>(globalData.changeOfBasisMatrices(2))),
        init::rDivM::view<3>::create(const_cast<real*>(globalData.changeOfBasisMatrices(3))),
      });
      for (int face = 0; face < 4; ++face) {
        auto& view = changeOfBasisMatricesView[face] ;
        for (int i = 0; i < view.shape(0); ++i) {
          for (int j = 0; j < view.shape(1); ++j) {
            view.shape(0);
            if (view.isInRange(i, j)) {
              // std::cout << "i,j = " << i << "\t" << j <<  "\t" << view.isInRange(i,j) << std::endl;
              view(i, j) *= filter->getFilterCoeff(j); // or i?!
            }
          }
        }
      }

      // volumeSum += self.db.kDivM[i][self.t('kl')] * self.I['lq'] * self.starMatrix(i)['qp']
      // Ignored: volumeSum += self.I['kq'] * self.sourceMatrix()['qp']

      auto stiffnessMatricesView = std::array<yateto::DenseTensorView<2, real>, 3>{
          init::kDivM::view<0>::create(const_cast<real*>(globalData.stiffnessMatrices(0))),
          init::kDivM::view<1>::create(const_cast<real*>(globalData.stiffnessMatrices(1))),
          init::kDivM::view<2>::create(const_cast<real*>(globalData.stiffnessMatrices(2))),
      };
      for (int d = 0; d < 3; ++d) {
        auto& view = stiffnessMatricesView[d] ;
        for (int i = 0; i < view.shape(0); ++i) {
          for (int j = 0; j < view.shape(1); ++j) {
            view.shape(0);
            if (view.isInRange(i, j)) {
              // std::cout << "i,j = " << i << "\t" << j <<  "\t" << view.isInRange(i,j) << std::endl;
              view(i, j) *= filter->getFilterCoeff(j); // or i?!
            }
          }
        }
      }
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

    void OnDevice::negateStiffnessMatrix(GlobalData &globalData) {
#ifdef ACL_DEVICE
      device::DeviceInstance& device = device::DeviceInstance::getInstance();
      for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
        const real scaleFactor = -1.0;
        device.algorithms.scaleArray(const_cast<real*>(globalData.stiffnessMatricesTransposed(transposedStiffness)),
                                     scaleFactor,
                                     init::kDivMT::size(transposedStiffness),
                                     device.api->getDefaultStream());
      }
#endif // ACL_DEVICE
    }
    void OnDevice::initSpecificGlobalData(GlobalData& globalData,
                                          memory::ManagedAllocator& allocator,
                                          CopyManagerT& copyManager,
                                          size_t alignment,
                                          seissol::memory::Memkind memkind) {
#ifdef ACL_DEVICE
      const size_t size = yateto::alignedUpper(tensor::replicateInitialLoadingM::size(),
                                               yateto::alignedReals<real>(alignment));
      real* plasticityStressReplication =
          static_cast<real*>(allocator.allocateMemory(size * sizeof(real),
                                                      alignment,
                                                      memkind));

      copyManager.template copyTensorToMemAndSetPtr<init::replicateInitialLoadingM>(plasticityStressReplication,
                                                                                    globalData.replicateStresses,
                                                                                    alignment);
#endif // ACL_DEVICE
    }

    real* OnDevice::DeviceCopyPolicy::copy(real const* first, real const* last, real*& mem) {
#ifdef ACL_DEVICE
      device::DeviceInstance& device = device::DeviceInstance::getInstance();
      const unsigned bytes = (last - first) * sizeof(real);
      device.api->copyTo(mem, first, bytes);
      mem += (last - first);
      return mem;
#else // ACL_DEVICE
      return nullptr;
#endif
    }

    void OnDevice::applyFilterMatrixToGlobalData(GlobalData& globalData, const seissol::kernels::Filter* filter) {
#ifdef ACL_DEVICE
      //TODO Implement on GPU
      // Be careful, need to consider pre-multiplied matrices, e.g. plusFluxMatrixAccessor!
#endif // ACL_DEVICE
    }
  } // namespace matrixmanip


template<typename MatrixManipPolicyT>
void GlobalDataInitializer<MatrixManipPolicyT>::init(GlobalData& globalData,
                                                     memory::ManagedAllocator& memoryAllocator,
                                                     enum seissol::memory::Memkind memkind) {
  MemoryProperties prop = MatrixManipPolicyT::getProperties();

  // We ensure that global matrices always start at an aligned memory address,
  // such that mixed cases with aligned and non-aligned global matrices do also work.
  unsigned globalMatrixMemSize = 0;
  globalMatrixMemSize += yateto::computeFamilySize<init::kDivM>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::kDivMT>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::rDivM>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::rT>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::fMrT>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::fP>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<nodal::init::V3mTo2nFace>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::project2nFaceTo3m>(yateto::alignedReals<real>(prop.alignment));

  globalMatrixMemSize += yateto::alignedUpper(tensor::evalAtQP::size(),  yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::alignedUpper(tensor::projectQP::size(), yateto::alignedReals<real>(prop.alignment));

#if defined(ACL_DEVICE) && defined(USE_PREMULTIPLY_FLUX)
  globalMatrixMemSize += yateto::computeFamilySize<init::plusFluxMatrices>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::minusFluxMatrices>(yateto::alignedReals<real>(prop.alignment));
#endif // ACL_DEVICE

  real* globalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory(globalMatrixMemSize * sizeof(real),
                                                                            prop.pagesizeHeap,
                                                                            memkind));

  real* globalMatrixMemPtr = globalMatrixMem;
  typename MatrixManipPolicyT::CopyManagerT copyManager;
  copyManager.template copyFamilyToMemAndSetPtr<init::kDivMT>(globalMatrixMemPtr, globalData.stiffnessMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::kDivM>(globalMatrixMemPtr, globalData.stiffnessMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::rDivM>(globalMatrixMemPtr, globalData.changeOfBasisMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::rT>(globalMatrixMemPtr, globalData.neighbourChangeOfBasisMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::fMrT>(globalMatrixMemPtr, globalData.localChangeOfBasisMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::fP>(globalMatrixMemPtr, globalData.neighbourFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<nodal::init::V3mTo2nFace>(globalMatrixMemPtr, globalData.V3mTo2nFace, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::project2nFaceTo3m>(globalMatrixMemPtr, globalData.project2nFaceTo3m, prop.alignment);

  copyManager.template copyTensorToMemAndSetPtr<init::evalAtQP>(globalMatrixMemPtr, globalData.evalAtQPMatrix, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<init::projectQP>(globalMatrixMemPtr, globalData.projectQPMatrix, prop.alignment);

#if defined(ACL_DEVICE) && defined(USE_PREMULTIPLY_FLUX)
  copyManager.template copyFamilyToMemAndSetPtr<init::plusFluxMatrices>(globalMatrixMemPtr, globalData.plusFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::minusFluxMatrices>(globalMatrixMemPtr, globalData.minusFluxMatrices, prop.alignment);
#endif // ACL_DEVICE

  assert(globalMatrixMemPtr == globalMatrixMem + globalMatrixMemSize);

  // @TODO Integrate this step into the code generator
  MatrixManipPolicyT::negateStiffnessMatrix(globalData);

  // Dynamic Rupture global matrices
  unsigned drGlobalMatrixMemSize = 0;
  drGlobalMatrixMemSize += yateto::computeFamilySize<init::V3mTo2nTWDivM>(yateto::alignedReals<real>(prop.alignment));
  drGlobalMatrixMemSize += yateto::computeFamilySize<init::V3mTo2n>(yateto::alignedReals<real>(prop.alignment));

  real* drGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory(drGlobalMatrixMemSize  * sizeof(real),
                                                                              prop.pagesizeHeap,
                                                                              memkind));

  real* drGlobalMatrixMemPtr = drGlobalMatrixMem;
  copyManager.template copyFamilyToMemAndSetPtr<init::V3mTo2nTWDivM>(drGlobalMatrixMemPtr, globalData.nodalFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::V3mTo2n>(drGlobalMatrixMemPtr, globalData.faceToNodalMatrices, prop.alignment);

  assert(drGlobalMatrixMemPtr == drGlobalMatrixMem + drGlobalMatrixMemSize);

  // Plasticity global matrices
  unsigned plasticityGlobalMatrixMemSize = 0;
  plasticityGlobalMatrixMemSize += yateto::alignedUpper(tensor::v::size(),    yateto::alignedReals<real>(prop.alignment));
  plasticityGlobalMatrixMemSize += yateto::alignedUpper(tensor::vInv::size(), yateto::alignedReals<real>(prop.alignment));

  real* plasticityGlobalMatrixMem
      = static_cast<real*>(memoryAllocator.allocateMemory(plasticityGlobalMatrixMemSize * sizeof(real),
                                                          prop.pagesizeHeap,
                                                          memkind));

  real* plasticityGlobalMatrixMemPtr = plasticityGlobalMatrixMem;
  copyManager.template copyTensorToMemAndSetPtr<init::v>(plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrix, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<init::vInv>(plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrixInverse, prop.alignment);

  assert(plasticityGlobalMatrixMemPtr == plasticityGlobalMatrixMem + plasticityGlobalMatrixMemSize);

  MatrixManipPolicyT::initSpecificGlobalData(globalData,
                                             memoryAllocator,
                                             copyManager,
                                             prop.pagesizeStack,
                                             memkind);


  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  auto filter = kernels::makeFilter(seissolParams.filter, 3);

  // TODO Implement for device
  MatrixManipPolicyT::applyFilterMatrixToGlobalData(globalData, filter.get());
}

template void GlobalDataInitializer<matrixmanip::OnHost>::init(GlobalData& globalData,
                                                               memory::ManagedAllocator& memoryAllocator,
                                                               enum memory::Memkind memkind);

template void GlobalDataInitializer<matrixmanip::OnDevice>::init(GlobalData& globalData,
                                                                 memory::ManagedAllocator& memoryAllocator,
                                                                 enum memory::Memkind memkind);


} // namespace seissol::initializers
