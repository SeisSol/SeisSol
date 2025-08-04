// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "GlobalData.h"
#include "GeneratedCode/init.h"
#include "Parallel/OpenMP.h"
#include <Common/ConfigHelper.h>
#include <DynamicRupture/FrictionLaws/TPCommon.h>
#include <DynamicRupture/Misc.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Memory/MemoryAllocator.h>
#include <cassert>
#include <cstddef>
#include <yateto.h>

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
template <typename Cfg>
struct OnHost {
  using CopyManagerT = typename yateto::DefaultCopyManager<real>;
  static MemoryProperties getProperties();
  static void negateStiffnessMatrix(GlobalDataCfg<Cfg>& globalData);
  static void initSpecificGlobalData(GlobalDataCfg<Cfg>& globalData,
                                     memory::ManagedAllocator& allocator,
                                     CopyManagerT& copyManager,
                                     size_t alignment,
                                     seissol::memory::Memkind memkind);
};

template <typename Cfg>
struct OnDevice {
  struct DeviceCopyPolicy {
    static real* copy(const real* first, const real* last, real*& mem);
  };
  using CopyManagerT = typename yateto::CopyManager<real, DeviceCopyPolicy>;
  static MemoryProperties getProperties();
  static void negateStiffnessMatrix(GlobalDataCfg<Cfg>& globalData);
  static void initSpecificGlobalData(GlobalDataCfg<Cfg>& globalData,
                                     memory::ManagedAllocator& allocator,
                                     CopyManagerT& copyManager,
                                     size_t alignment,
                                     seissol::memory::Memkind memkind);
};
} // namespace matrixmanip

// Generalized Global data initializers of SeisSol.
template <typename Cfg, typename MatrixManipPolicyT>
struct GlobalDataInitializer {
  static void init(GlobalDataCfg<Cfg>& globalData,
                   memory::ManagedAllocator& memoryAllocator,
                   enum memory::Memkind memkind);
};
} // namespace seissol::initializer

namespace seissol::initializer {
namespace matrixmanip {
template <typename Cfg>
MemoryProperties OnHost<Cfg>::getProperties() {
  // returns MemoryProperties initialized with default values i.e., CPU memory properties
  return {};
}

template <typename Cfg>
void OnHost<Cfg>::negateStiffnessMatrix(GlobalDataCfg<Cfg>& globalData) {
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    real* matrix = const_cast<real*>(globalData.stiffnessMatricesTransposed(transposedStiffness));
    for (unsigned i = 0; i < init::kDivMT<Cfg>::size(transposedStiffness); ++i) {
      matrix[i] *= -1.0;
    }
  }
}

template <typename Cfg>
void OnHost<Cfg>::initSpecificGlobalData(GlobalDataCfg<Cfg>& globalData,
                                         memory::ManagedAllocator& allocator,
                                         CopyManagerT& copyManager,
                                         size_t alignment,
                                         seissol::memory::Memkind memkind) {
  // thread-local LTS integration buffers
  const auto numThreads = OpenMP::threadCount();
  const auto allocSize = 4 * static_cast<std::size_t>(tensor::I<Cfg>::size());
  auto* integrationBufferLTS = reinterpret_cast<real*>(
      allocator.allocateMemory(numThreads * allocSize * sizeof(real), alignment, memkind));

// initialize w.r.t. NUMA
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    const auto threadOffset = OpenMP::threadId() * allocSize;
    for (std::size_t dof = 0; dof < allocSize; ++dof) {
      integrationBufferLTS[dof + threadOffset] = static_cast<real>(0.0);
    }
  }

  globalData.integrationBufferLTS = integrationBufferLTS;
}

template <typename Cfg>
MemoryProperties OnDevice<Cfg>::getProperties() {
  MemoryProperties prop{};
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  prop.alignment = device.api->getGlobMemAlignment();
  prop.pagesizeHeap = prop.alignment;
  prop.pagesizeStack = prop.alignment;
#endif
  return prop;
}

template <typename Cfg>
void OnDevice<Cfg>::negateStiffnessMatrix(GlobalDataCfg<Cfg>& globalData) {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    const real scaleFactor = -1.0;
    device.algorithms.scaleArray(
        const_cast<real*>(globalData.stiffnessMatricesTransposed(transposedStiffness)),
        scaleFactor,
        init::kDivMT<Cfg>::size(transposedStiffness),
        device.api->getDefaultStream());
  }
#endif // ACL_DEVICE
}

template <typename Cfg>
void OnDevice<Cfg>::initSpecificGlobalData(GlobalDataCfg<Cfg>& globalData,
                                           memory::ManagedAllocator& allocator,
                                           CopyManagerT& copyManager,
                                           size_t alignment,
                                           seissol::memory::Memkind memkind) {
#ifdef ACL_DEVICE
  const size_t size = yateto::alignedUpper(tensor::replicateInitialLoadingM<Cfg>::size(),
                                           yateto::alignedReals<real>(alignment));
  real* plasticityStressReplication =
      static_cast<real*>(allocator.allocateMemory(size * sizeof(real), alignment, memkind));

  copyManager.template copyTensorToMemAndSetPtr<init::replicateInitialLoadingM<Cfg>>(
      plasticityStressReplication, globalData.replicateStresses, alignment);
#endif // ACL_DEVICE
}

template <typename Cfg>
real* OnDevice<Cfg>::DeviceCopyPolicy::copy(const real* first, const real* last, real*& mem) {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  const std::size_t bytes = (last - first) * sizeof(real);
  device.api->copyTo(mem, first, bytes);
  mem += (last - first);
  return mem;
#else // ACL_DEVICE
  return nullptr;
#endif
}

} // namespace matrixmanip

template <typename Cfg, typename MatrixManipPolicyT>
void GlobalDataInitializer<Cfg, MatrixManipPolicyT>::init(GlobalDataCfg<Cfg>& globalData,
                                                          memory::ManagedAllocator& memoryAllocator,
                                                          enum seissol::memory::Memkind memkind) {
  const MemoryProperties prop = MatrixManipPolicyT::getProperties();

  // We ensure that global matrices always start at an aligned memory address,
  // such that mixed cases with aligned and non-aligned global matrices do also work.
  unsigned globalMatrixMemSize = 0;
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::kDivM<Cfg>>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::kDivMT<Cfg>>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::rDivM<Cfg>>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::rT<Cfg>>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::fMrT<Cfg>>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::computeFamilySize<init::fP<Cfg>>(yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<nodal::init::V3mTo2nFace<Cfg>>(
      yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::project2nFaceTo3m<Cfg>>(
      yateto::alignedReals<real>(prop.alignment));

  globalMatrixMemSize += yateto::alignedUpper(tensor::evalAtQP<Cfg>::size(),
                                              yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::alignedUpper(tensor::projectQP<Cfg>::size(),
                                              yateto::alignedReals<real>(prop.alignment));

#ifdef USE_VISCOELASTIC2
  globalMatrixMemSize += yateto::alignedUpper(tensor::selectAne<Cfg>::size(),
                                              yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::alignedUpper(tensor::selectEla<Cfg>::size(),
                                              yateto::alignedReals<real>(prop.alignment));
#endif

#if defined(ACL_DEVICE) && defined(USE_PREMULTIPLY_FLUX)
  globalMatrixMemSize += yateto::computeFamilySize<init::plusFluxMatrices<Cfg>>(
      yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::computeFamilySize<init::minusFluxMatrices<Cfg>>(
      yateto::alignedReals<real>(prop.alignment));
#endif // ACL_DEVICE

  globalMatrixMemSize += yateto::alignedUpper(tensor::resample<Cfg>::size(),
                                              yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize += yateto::alignedUpper(tensor::quadweights<Cfg>::size(),
                                              yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::alignedUpper(dr::misc::NumTpGridPoints, yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::alignedUpper(dr::misc::NumTpGridPoints, yateto::alignedReals<real>(prop.alignment));
  globalMatrixMemSize +=
      yateto::alignedUpper(dr::misc::NumTpGridPoints, yateto::alignedReals<real>(prop.alignment));

  real* globalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory(
      globalMatrixMemSize * sizeof(real), prop.pagesizeHeap, memkind));

  real* globalMatrixMemPtr = globalMatrixMem;
  typename MatrixManipPolicyT::CopyManagerT copyManager;
  copyManager.template copyFamilyToMemAndSetPtr<init::kDivMT<Cfg>>(
      globalMatrixMemPtr, globalData.stiffnessMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::kDivM<Cfg>>(
      globalMatrixMemPtr, globalData.stiffnessMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::rDivM<Cfg>>(
      globalMatrixMemPtr, globalData.changeOfBasisMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::rT<Cfg>>(
      globalMatrixMemPtr, globalData.neighborChangeOfBasisMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::fMrT<Cfg>>(
      globalMatrixMemPtr, globalData.localChangeOfBasisMatricesTransposed, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::fP<Cfg>>(
      globalMatrixMemPtr, globalData.neighborFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<nodal::init::V3mTo2nFace<Cfg>>(
      globalMatrixMemPtr, globalData.v3mTo2nFace, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::project2nFaceTo3m<Cfg>>(
      globalMatrixMemPtr, globalData.project2nFaceTo3m, prop.alignment);

  copyManager.template copyTensorToMemAndSetPtr<init::evalAtQP<Cfg>>(
      globalMatrixMemPtr, globalData.evalAtQPMatrix, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<init::projectQP<Cfg>>(
      globalMatrixMemPtr, globalData.projectQPMatrix, prop.alignment);

#ifdef USE_VISCOELASTIC2
  copyManager.template copyTensorToMemAndSetPtr<init::selectAne<Cfg>>(
      globalMatrixMemPtr, globalData.selectAne, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<init::selectEla<Cfg>>(
      globalMatrixMemPtr, globalData.selectEla, prop.alignment);
#endif

#if defined(ACL_DEVICE) && defined(USE_PREMULTIPLY_FLUX)
  copyManager.template copyFamilyToMemAndSetPtr<init::plusFluxMatrices<Cfg>>(
      globalMatrixMemPtr, globalData.plusFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::minusFluxMatrices<Cfg>>(
      globalMatrixMemPtr, globalData.minusFluxMatrices, prop.alignment);
#endif // ACL_DEVICE

  copyManager.template copyTensorToMemAndSetPtr<init::resample<Cfg>>(
      globalMatrixMemPtr, globalData.resampleMatrix, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<init::quadweights<Cfg>>(
      globalMatrixMemPtr, globalData.spaceWeights, prop.alignment);

  // a bit more manual
  {
    const auto data =
        seissol::dr::friction_law::tp::InverseFourierCoefficients<dr::misc::NumTpGridPoints>();
    globalData.tpInverseFourierCoefficients = globalMatrixMemPtr;
    globalMatrixMemPtr +=
        yateto::alignedUpper(dr::misc::NumTpGridPoints, yateto::alignedReals<real>(prop.alignment));
    seissol::memory::memcopyTyped<real>(globalData.tpInverseFourierCoefficients,
                                        data.data().data(),
                                        dr::misc::NumTpGridPoints,
                                        memkind,
                                        memory::Standard);
  }

  {
    const auto data = seissol::dr::friction_law::tp::GridPoints<dr::misc::NumTpGridPoints>();
    globalData.tpGridPoints = globalMatrixMemPtr;
    globalMatrixMemPtr +=
        yateto::alignedUpper(dr::misc::NumTpGridPoints, yateto::alignedReals<real>(prop.alignment));
    seissol::memory::memcopyTyped<real>(globalData.tpGridPoints,
                                        data.data().data(),
                                        dr::misc::NumTpGridPoints,
                                        memkind,
                                        memory::Standard);
  }

  {
    const auto data =
        seissol::dr::friction_law::tp::GaussianHeatSource<dr::misc::NumTpGridPoints>();
    globalData.heatSource = globalMatrixMemPtr;
    globalMatrixMemPtr +=
        yateto::alignedUpper(dr::misc::NumTpGridPoints, yateto::alignedReals<real>(prop.alignment));

    seissol::memory::memcopyTyped<real>(globalData.heatSource,
                                        data.data().data(),
                                        dr::misc::NumTpGridPoints,
                                        memkind,
                                        memory::Standard);
  }

  assert(globalMatrixMemPtr == globalMatrixMem + globalMatrixMemSize);

  // @TODO Integrate this step into the code generator
  MatrixManipPolicyT::negateStiffnessMatrix(globalData);

  // Dynamic Rupture global matrices
  unsigned drGlobalMatrixMemSize = 0;
  drGlobalMatrixMemSize += yateto::computeFamilySize<init::V3mTo2nTWDivM<Cfg>>(
      yateto::alignedReals<real>(prop.alignment));
  drGlobalMatrixMemSize +=
      yateto::computeFamilySize<init::V3mTo2n<Cfg>>(yateto::alignedReals<real>(prop.alignment));

  real* drGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory(
      drGlobalMatrixMemSize * sizeof(real), prop.pagesizeHeap, memkind));

  real* drGlobalMatrixMemPtr = drGlobalMatrixMem;
  copyManager.template copyFamilyToMemAndSetPtr<init::V3mTo2nTWDivM<Cfg>>(
      drGlobalMatrixMemPtr, globalData.nodalFluxMatrices, prop.alignment);
  copyManager.template copyFamilyToMemAndSetPtr<init::V3mTo2n<Cfg>>(
      drGlobalMatrixMemPtr, globalData.faceToNodalMatrices, prop.alignment);

  assert(drGlobalMatrixMemPtr == drGlobalMatrixMem + drGlobalMatrixMemSize);

  // Plasticity global matrices
  unsigned plasticityGlobalMatrixMemSize = 0;
  plasticityGlobalMatrixMemSize +=
      yateto::alignedUpper(tensor::v<Cfg>::size(), yateto::alignedReals<real>(prop.alignment));
  plasticityGlobalMatrixMemSize +=
      yateto::alignedUpper(tensor::vInv<Cfg>::size(), yateto::alignedReals<real>(prop.alignment));

  real* plasticityGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory(
      plasticityGlobalMatrixMemSize * sizeof(real), prop.pagesizeHeap, memkind));

  real* plasticityGlobalMatrixMemPtr = plasticityGlobalMatrixMem;
  copyManager.template copyTensorToMemAndSetPtr<init::v<Cfg>>(
      plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrix, prop.alignment);
  copyManager.template copyTensorToMemAndSetPtr<init::vInv<Cfg>>(
      plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrixInverse, prop.alignment);

  assert(plasticityGlobalMatrixMemPtr == plasticityGlobalMatrixMem + plasticityGlobalMatrixMemSize);

  MatrixManipPolicyT::initSpecificGlobalData(
      globalData, memoryAllocator, copyManager, prop.pagesizeStack, memkind);
}

} // namespace seissol::initializer

namespace seissol {

void GlobalData::init(std::size_t configId) {
  std::visit(
      [&](auto cfg) {
        using Cfg = decltype(cfg);

        GlobalDataCfg<Cfg> host{};
        initializer::GlobalDataInitializer<Cfg, initializer::matrixmanip::OnHost<Cfg>>::init(
            host, allocator, memkindHost);
        std::get<std::optional<GlobalDataCfg<Cfg>>>(onHost) = host;

        if constexpr (isDeviceOn()) {
          GlobalDataCfg<Cfg> device{};
          initializer::GlobalDataInitializer<Cfg, initializer::matrixmanip::OnDevice<Cfg>>::init(
              device, allocator, memkindDevice);
          std::get<std::optional<GlobalDataCfg<Cfg>>>(onDevice) = device;
        }
      },
      ConfigVariantList[configId]);
}

} // namespace seissol
