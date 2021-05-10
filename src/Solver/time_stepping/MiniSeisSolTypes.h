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

struct dofs {
  using type = std::array<real, tensor::Q::size()>;
};


// TODO(Lukas) Find good naming scheme for structs
// maybe lowercase?
struct buffer {
  using type = real*;
};

struct derivatives {
  using type = real*;
};

struct cellLocalInformation {
  using type = CellLocalInformation;
};

struct faceNeighbors {
  using type = std::array<real*, 4>;
};

// TODO(Lukas) Ugh... name
struct localIntegrationData {
  using type = LocalIntegrationData;
};

struct neighborIntegrationData {
  using type = NeighboringIntegrationData;
};

struct material {
  using type = CellMaterialData;
};


struct plasticity {
  using type = PlasticityData;
};

struct cellDrMapping {
  using type = CellDRMapping;
};

struct BoundaryMapping {
  using type = CellBoundaryMapping;
};

struct Pstrain {
  using type = real[7];
};

struct displacements {
  using type = real*;
};


struct BuffersBucket {
  using type = real;
};

struct DisplacementsBucket {
  using type = real;
};

using element_storage_t = mneme::MultiStorage<
    mneme::DataLayout::SoA,
    dofs,
    buffer,
    derivatives,
    cellLocalInformation,
    faceNeighbors,
    localIntegrationData,
    neighborIntegrationData,
    displacements
>;

// Buckets
using buffers_bucket_storage_t = mneme::SingleStorage<BuffersBucket>;
using buffers_bucket_displacements_t = mneme::SingleStorage<DisplacementsBucket>;

// TODO(Lukas) Temporary
struct ProxyData {
  ProxyData(
      std::shared_ptr<element_storage_t> elementStorage,
      mneme::LayeredPlan<seissol::InteriorLayer, seissol::CopyLayer, seissol::GhostLayer> elementStoragePlan,
      std::shared_ptr<buffers_bucket_storage_t> buffersBucket,
      mneme::LayeredPlan<seissol::InteriorLayer, seissol::CopyLayer, seissol::GhostLayer> buffersBucketPlan
      )
      :
      elementStorage(std::move(elementStorage)),
      buffersBucket(std::move(buffersBucket)),
      elementStoragePlan(std::move(elementStoragePlan)),
      buffersBucketPlan(std::move(buffersBucketPlan))
      {

  }

  std::shared_ptr<element_storage_t> elementStorage;
  mneme::LayeredPlan<seissol::InteriorLayer, seissol::CopyLayer, seissol::GhostLayer> elementStoragePlan;

  std::shared_ptr<buffers_bucket_storage_t> buffersBucket;
  mneme::LayeredPlan<seissol::InteriorLayer, seissol::CopyLayer, seissol::GhostLayer> buffersBucketPlan;
  //std::shared_ptr<buffers_bucket_displacements_t> displacementBucket;

};

}

#endif //SEISSOL_MINISEISSOLTYPES_H
