#ifndef SEISSOL_MINISEISSOLTYPES_H
#define SEISSOL_MINISEISSOLTYPES_H


#include "Initializer/MemoryManager.h"

#include "Kernels/precision.hpp"
#include <mneme/storage.hpp>
#include "mneme/plan.hpp"
#include "mneme/view.hpp"
// TODO(Lukas) Names!
namespace seissol {
struct GhostLayer : public mneme::Layer {
};
struct InteriorLayer : public mneme::Layer {
};
struct CopyLayer : public mneme::Layer {
};

template<typename type>
using standardAllocator = mneme::AlignedAllocator<type, ALIGNMENT>;

struct dofs {
  using type = std::array<real, tensor::Q::size()>;
  using allocator = standardAllocator<type>;
};


// TODO(Lukas) Find good naming scheme for structs
// maybe lowercase?
struct buffer {
  using type = real*;
  using allocator = standardAllocator<type>;
};

struct derivatives {
  using type = real*;
  using allocator = standardAllocator<type>;
};

struct cellLocalInformation {
  using type = CellLocalInformation;
  using allocator = standardAllocator<type>;
};

struct faceNeighbors {
  using type = std::array<real*, 4>;
  using allocator = standardAllocator<type>;
};

// TODO(Lukas) Ugh... name
struct localIntegrationData {
  using type = LocalIntegrationData;
  using allocator = standardAllocator<type>;
};

struct neighborIntegrationData {
  using type = NeighboringIntegrationData;
  using allocator = standardAllocator<type>;
};

struct material {
  using type = CellMaterialData;
  using allocator = standardAllocator<type>;
};


struct plasticity {
  using type = PlasticityData;
  using allocator = standardAllocator<type>;
};

struct cellDrMapping {
  using type = std::array<CellDRMapping, 4>;
  using allocator = standardAllocator<type>;
};

struct BoundaryMapping {
  using type = CellBoundaryMapping;
  using allocator = standardAllocator<type>;
};

/*struct Pstrain {
  using type = std::array<real,7>;
  using allocator = standardAllocator<type>;
};
 */

struct displacements {
  using type = real*;
  using allocator = standardAllocator<type>;
};


struct BuffersBucket {
  using type = real;
  using allocator = standardAllocator<type>;
};

struct DisplacementsBucket {
  using type = real;
  using allocator = standardAllocator<type>;
};

using element_storage_t = mneme::MultiStorage<
    mneme::DataLayout::SoA,
    dofs,
    buffer,
    derivatives,
    cellLocalInformation,
    cellDrMapping,
    faceNeighbors,
    localIntegrationData,
    neighborIntegrationData,
    displacements
>;

// Buckets
using buffers_bucket_storage_t = mneme::SingleStorage<BuffersBucket>;
using buffers_bucket_displacements_t = mneme::SingleStorage<DisplacementsBucket>;

struct timeDerivativePlus {
  using type = real*;
  using allocator = standardAllocator<type>;
};
struct timeDerivativeMinus {
  using type = real*;
  using allocator = standardAllocator<type>;
};

struct imposedStatePlus {
  using type = std::array<real, tensor::QInterpolated::size()>;
  using allocator = standardAllocator<type>;
};

struct imposedStateMinus {
  using type = std::array<real, tensor::QInterpolated::size()>;
  using allocator = standardAllocator<type>;
};

struct godunovData {
  using type = DRGodunovData;
  using allocator = standardAllocator<type>;
};

struct fluxSolverPlus {
  using type = std::array<real, tensor::fluxSolver::size()>;
  using allocator = standardAllocator<type>;
};

struct fluxSolverMinus {
  using type = std::array<real, tensor::fluxSolver::size()>;
  using allocator = standardAllocator<type>;
};

struct faceInformation {
  using type = DRFaceInformation;
  using allocator = standardAllocator<type>;
};

struct waveSpeedsPlus {
  using type = model::IsotropicWaveSpeeds;
  using allocator = standardAllocator<type>;
};

struct waveSpeedsMinus {
  using type = model::IsotropicWaveSpeeds;
  using allocator = standardAllocator<type>;
};

using dynamic_rupture_storage_t = mneme::MultiStorage<
    mneme::DataLayout::SoA,
    timeDerivativePlus,
    timeDerivativeMinus,
    imposedStatePlus,
    imposedStateMinus,
    godunovData,
    fluxSolverPlus,
    fluxSolverMinus,
    faceInformation,
    waveSpeedsPlus,
    waveSpeedsMinus
    >;

// TODO(Lukas) Temporary
struct ProxyData {
  using plan_t = mneme::CombinedLayeredPlan<seissol::GhostLayer, seissol::CopyLayer, seissol::InteriorLayer>;

  ProxyData(
      std::shared_ptr<element_storage_t> elementStorage,
      plan_t elementStoragePlan,
      std::shared_ptr<buffers_bucket_storage_t> buffersBucket,
      plan_t buffersBucketPlan,
      std::shared_ptr<dynamic_rupture_storage_t> dynamicRuptureStorage,
      plan_t dynamicRuptureStoragePlan
      )
      :
      elementStorage(std::move(elementStorage)),
      buffersBucket(std::move(buffersBucket)),
      elementStoragePlan(std::move(elementStoragePlan)),
      buffersBucketPlan(std::move(buffersBucketPlan)),
      dynamicRuptureStorage(std::move(dynamicRuptureStorage)),
      dynamicRuptureStoragePlan(std::move(dynamicRuptureStoragePlan))
      {


  }

  auto getElementView() {
    return mneme::createViewFactory()
        .withPlan(elementStoragePlan)
        .withStorage(elementStorage)
        .withClusterId(0)
        .createDenseView<InteriorLayer>();
  }

  auto getDynamicRuptureView() {
    return mneme::createViewFactory()
        .withPlan(dynamicRuptureStoragePlan)
        .withStorage(dynamicRuptureStorage)
        .withClusterId(0)
        .createDenseView<InteriorLayer>();
  }

  std::shared_ptr<element_storage_t> elementStorage;
  plan_t elementStoragePlan;

  std::shared_ptr<buffers_bucket_storage_t> buffersBucket;
  plan_t buffersBucketPlan;
  //std::shared_ptr<buffers_bucket_displacements_t> displacementBucket;

  std::shared_ptr<dynamic_rupture_storage_t> dynamicRuptureStorage;
  plan_t dynamicRuptureStoragePlan;
};

}

#endif //SEISSOL_MINISEISSOLTYPES_H
