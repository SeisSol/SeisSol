// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "GlobalData.h"
#include "generated_code/init.h"
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Memory/MemoryAllocator.h>
#include <cassert>
#include <cstddef>
#include <tensor.h>
#include <yateto.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace seissol::initializer {
namespace matrixmanip {
MemoryProperties OnHost::getProperties() {
  // returns MemoryProperties initialized with default values i.e., CPU memory properties
  return {};
}

void OnHost::negateStiffnessMatrix(GlobalData& globalData) {
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    real* matrix = const_cast<real*>(globalData.stiffnessMatricesTransposed(transposedStiffness));
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
  int numberOfThreads = 1;
#ifdef _OPENMP
  numberOfThreads = omp_get_max_threads();
#endif
  auto* integrationBufferLTS = reinterpret_cast<real*>(allocator.allocateMemory(
      numberOfThreads * (4 * tensor::I::size()) * sizeof(real), alignment, memkind));

// initialize w.r.t. NUMA
#ifdef _OPENMP
#pragma omp parallel
  {
    const std::size_t threadOffset = omp_get_thread_num() * (4 * tensor::I::size());
#else
  std::size_t threadOffset = 0;
#endif
    for (std::size_t dof = 0; dof < (4 * tensor::I::size()); ++dof) {
      integrationBufferLTS[dof + threadOffset] = (real)0.0;
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

void OnDevice::negateStiffnessMatrix(GlobalData& globalData) {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    const real scaleFactor = -1.0;
    device.algorithms.scaleArray(
        const_cast<real*>(globalData.stiffnessMatricesTransposed(transposedStiffness)),
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
      static_cast<real*>(allocator.allocateMemory(size * sizeof(real), alignment, memkind));

  copyManager.template copyTensorToMemAndSetPtr<init::replicateInitialLoadingM>(
      plasticityStressReplication, globalData.replicateStresses, alignment);
#endif // ACL_DEVICE
}

real* OnDevice::DeviceCopyPolicy::copy(const real* first, const real* last, real*& mem) {
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

} // namespace matrixmanip

template <typename MatrixManipPolicyT>
void GlobalDataInitializer<MatrixManipPolicyT>::init(GlobalData& globalData,
                                                     memory::ManagedAllocator& memoryAllocator,
                                                     enum seissol::memory::Memkind memkind) {
  const MemoryProperties prop = MatrixManipPolicyT::getProperties();

  // We ensure that global matrices always start at an aligned memory address,
  // such that mixed cases with aligned and non-aligned global matrices do also work.
  unsigned globalMatrixMemSize = 0;
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::kDivM>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::kDivMT>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::rDivM>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::rT>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::fMrT>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::fP>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<nodal::init::V3mTo2nFace>(
      yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::project2nFaceTo3m>(
      yateto::alignedReals<real>(prop.alignment));

  globalMatrixMemSize +=
      yateto::alignedUpper(tensor::evalAtQP::size(), yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::alignedUpper(tensor::projectQP::size(), yateto::alignedReals<real>(prop.alignment));

#if defined(ACL_DEVICE) && defined(USE_PREMULTIPLY_FLUX)
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::plusFluxMatrices>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::minusFluxMatrices>(
      yateto::alignedReals<real>(prop.alignment));
#endif // ACL_DEVICE

  real* globalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory(
      globalMatrixMemSize * sizeof(real), prop.pagesizeHeap, memkind));

  real* globalMatrixMemPtr = globalMatrixMem;
  typename MatrixManipPolicyT::CopyManagerT copyManager;
  copyManager.template copyFamilyToMemAndSetPtr<init::kDivMT>(
      globalMatrixMemPtr, globalData.stiffnessMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::kDivM>(
      globalMatrixMemPtr, globalData.stiffnessMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::rDivM>(
      globalMatrixMemPtr, globalData.changeOfBasisMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::rT>(
      globalMatrixMemPtr, globalData.neighbourChangeOfBasisMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::fMrT>(
      globalMatrixMemPtr, globalData.localChangeOfBasisMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::fP>(
      globalMatrixMemPtr, globalData.neighbourFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<nodal::init::V3mTo2nFace>(
      globalMatrixMemPtr, globalData.V3mTo2nFace, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::project2nFaceTo3m>(
      globalMatrixMemPtr, globalData.project2nFaceTo3m, prop.alignment);

  copyManager.template copyTensorToMemAndSetPtr<init::evalAtQP>(
      globalMatrixMemPtr, globalData.evalAtQPMatrix, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<init::projectQP>(
      globalMatrixMemPtr, globalData.projectQPMatrix, prop.alignment);

#if defined(ACL_DEVICE) && defined(USE_PREMULTIPLY_FLUX)
  copyManager.template copyFamilyToMemAndSetPtr<init::plusFluxMatrices>(
      globalMatrixMemPtr, globalData.plusFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::minusFluxMatrices>(
      globalMatrixMemPtr, globalData.minusFluxMatrices, prop.alignment);
#endif // ACL_DEVICE

  assert(globalMatrixMemPtr == globalMatrixMem + globalMatrixMemSize);

  // @TODO Integrate this step into the code generator
  MatrixManipPolicyT::negateStiffnessMatrix(globalData);

  // Dynamic Rupture global matrices
  unsigned drGlobalMatrixMemSize = 0;
  drGlobalMatrixMemSize +=
      yateto::computeFamilySize<init::V3mTo2nTWDivM>(yateto::alignedReals<real>(prop.alignment));
  drGlobalMatrixMemSize +=
      yateto::computeFamilySize<init::V3mTo2n>(yateto::alignedReals<real>(prop.alignment));

  real* drGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory(
      drGlobalMatrixMemSize * sizeof(real), prop.pagesizeHeap, memkind));

  real* drGlobalMatrixMemPtr = drGlobalMatrixMem;
  copyManager.template copyFamilyToMemAndSetPtr<init::V3mTo2nTWDivM>(
      drGlobalMatrixMemPtr, globalData.nodalFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::V3mTo2n>(
      drGlobalMatrixMemPtr, globalData.faceToNodalMatrices, prop.alignment);

  assert(drGlobalMatrixMemPtr == drGlobalMatrixMem + drGlobalMatrixMemSize);

  // Plasticity global matrices
  unsigned plasticityGlobalMatrixMemSize = 0;
  plasticityGlobalMatrixMemSize +=
      yateto::alignedUpper(tensor::v::size(), yateto::alignedReals<real>(prop.alignment));
  plasticityGlobalMatrixMemSize +=
      yateto::alignedUpper(tensor::vInv::size(), yateto::alignedReals<real>(prop.alignment));

  real* plasticityGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory(
      plasticityGlobalMatrixMemSize * sizeof(real), prop.pagesizeHeap, memkind));

  real* plasticityGlobalMatrixMemPtr = plasticityGlobalMatrixMem;
  copyManager.template copyTensorToMemAndSetPtr<init::v>(
      plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrix, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<init::vInv>(
      plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrixInverse, prop.alignment);

  assert(plasticityGlobalMatrixMemPtr == plasticityGlobalMatrixMem + plasticityGlobalMatrixMemSize);

  MatrixManipPolicyT::initSpecificGlobalData(
      globalData, memoryAllocator, copyManager, prop.pagesizeStack, memkind);
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
