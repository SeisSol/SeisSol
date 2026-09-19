// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.h"

#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Output/OutputAux.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Helper.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <yateto.h>

namespace seissol::dr::output {
void ReceiverBasedOutputBuilder::setMeshReader(const seissol::geometry::MeshReader* reader) {
  meshReader_ = reader;
  localRank_ = Mpi::mpi.rank();
}

void ReceiverBasedOutputBuilder::setLtsData(LTS::Storage& userWpStorage,
                                            LTS::Backmap& userWpBackmap,
                                            DynamicRupture::Storage& userDrStorage) {
  wpStorage_ = &userWpStorage;
  wpBackmap_ = &userWpBackmap;
  drStorage_ = &userDrStorage;
}

void ReceiverBasedOutputBuilder::setVariableList(const std::vector<std::size_t>& variables) {
  this->variables_ = variables;
}

void ReceiverBasedOutputBuilder::setFaceToLtsMap(
    ::seissol::initializer::StorageBackmap<1>* faceToLtsMap) {
  this->faceToLtsMap_ = faceToLtsMap;
}

namespace {
struct GhostElement {
  std::pair<std::size_t, int> data;
  std::size_t index{};
};

template <typename T1, typename T2>
struct HashPair {
  std::size_t operator()(const std::pair<T1, T2>& data) const {
    // Taken from: https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x
    // (probably any other lcg-like hash function would work as well)
    const std::hash<T1> hasher1;
    const std::hash<T2> hasher2;
    std::size_t seed = hasher1(data.first);
    seed ^= hasher2(data.second) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
    return seed;
  }
};
} // namespace

void ReceiverBasedOutputBuilder::initTopology() {
  auto& topology = outputData_->topology;

  // Bucket the receivers by face and, within a face, by position. Buckets are created on first
  // touch, so the order the generators produced is preserved: face-major for the elementwise
  // output, file order for the pickpoint output.
  struct PointBucket {
    std::vector<std::size_t> receiverIds;
    std::size_t nearestGpIndex{};
    std::size_t nearestInternalGpIndex{};
  };
  struct FaceBucket {
    std::size_t faultFaceIndex{};
    std::size_t elementIndex{};
    std::size_t localFaceSideId{};
    std::vector<PointBucket> points;
    std::map<std::array<double, 3>, std::size_t> pointIds;
  };

  std::vector<FaceBucket> faceBuckets;
  std::unordered_map<std::size_t, std::size_t> faceIds;

  for (const auto& [index, receiver] : common::enumerate(outputData_->receivers)) {
    if (!receiver.isInside) {
      continue;
    }

    const auto faceKey = faceToLtsMap_->get(receiver.faultFaceIndex).global;
    if (faceIds.find(faceKey) == faceIds.end()) {
      faceIds[faceKey] = faceBuckets.size();
      auto& bucket = faceBuckets.emplace_back();
      bucket.faultFaceIndex = static_cast<std::size_t>(receiver.faultFaceIndex);
      bucket.elementIndex = static_cast<std::size_t>(receiver.elementIndex);
      bucket.localFaceSideId = static_cast<std::size_t>(receiver.localFaceSideId);
    }
    auto& faceBucket = faceBuckets[faceIds.at(faceKey)];

    // Exact comparison is intended: receivers which share a position are generated from the very
    // same coordinates -- fused simulations, or a point given twice in the parameter file. The
    // lookup is scoped to the face, since the node points of two adjacent faces may coincide
    // bit-for-bit (a corner of the reference face maps to a mesh vertex exactly).
    const std::array<double, 3> coords{
        receiver.global.coords[0], receiver.global.coords[1], receiver.global.coords[2]};
    if (faceBucket.pointIds.find(coords) == faceBucket.pointIds.end()) {
      faceBucket.pointIds[coords] = faceBucket.points.size();
      auto& bucket = faceBucket.points.emplace_back();
      bucket.nearestGpIndex = static_cast<std::size_t>(receiver.nearestGpIndex);
      bucket.nearestInternalGpIndex = static_cast<std::size_t>(receiver.nearestInternalGpIndex);
    }
    faceBucket.points[faceBucket.pointIds.at(coords)].receiverIds.push_back(index);
  }

  // Flatten the buckets, renumbering the receivers on the way. This is what makes the point range
  // of a face and the receiver range of a point contiguous, so that the topology needs offsets
  // only.
  std::vector<Receiver> receivers;
  receivers.reserve(outputData_->receivers.size());

  for (const auto& faceBucket : faceBuckets) {
    OutputFace face{};
    face.faultFaceIndex = faceBucket.faultFaceIndex;
    face.position = faceToLtsMap_->get(faceBucket.faultFaceIndex);
    face.elementIndex = faceBucket.elementIndex;
    face.localFaceSideId = faceBucket.localFaceSideId;
    topology.addFace(face);

    for (const auto& pointBucket : faceBucket.points) {
      OutputPoint point{};
      point.faceId = topology.faceCount() - 1;
      point.nearestGpIndex = pointBucket.nearestGpIndex;
      point.nearestInternalGpIndex = pointBucket.nearestInternalGpIndex;
      topology.addPoint(point, pointBucket.receiverIds.size());

      for (const auto receiverId : pointBucket.receiverIds) {
        receivers.push_back(outputData_->receivers[receiverId]);
      }
    }
  }

  outputData_->receivers = std::move(receivers);
  assert(outputData_->receivers.size() == topology.receiverCount());
}

void ReceiverBasedOutputBuilder::initBasisFunctions() {
  const auto& faultInfo = meshReader_->getFault();
  const auto& elementsInfo = meshReader_->getElements();
  const auto& verticesInfo = meshReader_->getVertices();
  const auto& mpiGhostMetadata = meshReader_->getGhostlayerMetadata();

  auto& topology = outputData_->topology;

  for (std::size_t faceId = 0; faceId < topology.faceCount(); ++faceId) {
    const auto& face = topology.faces[faceId];
    const auto elementIndex = faultInfo[face.faultFaceIndex].element;
    const auto& element = elementsInfo[elementIndex];
    const auto neighborElementIndex = faultInfo[face.faultFaceIndex].neighborElement;

    const VrtxCoords* elemCoords[Cell::NumVertices]{};
    for (size_t vertexIdx = 0; vertexIdx < Cell::NumVertices; ++vertexIdx) {
      const auto address = element.vertices[vertexIdx];
      elemCoords[vertexIdx] = &(verticesInfo[address].coords);
    }

    const VrtxCoords* neighborElemCoords[Cell::NumVertices]{};
    if (neighborElementIndex >= 0) {
      for (size_t vertexIdx = 0; vertexIdx < Cell::NumVertices; ++vertexIdx) {
        const auto address = elementsInfo[neighborElementIndex].vertices[vertexIdx];
        neighborElemCoords[vertexIdx] = &(verticesInfo[address].coords);
      }
    } else {
      const auto faultSide = faultInfo[face.faultFaceIndex].side;
      const auto neighborRank = element.neighborRanks[faultSide];
      const auto& ghostMetadataItr = mpiGhostMetadata.find(neighborRank);
      assert(ghostMetadataItr != mpiGhostMetadata.end());

      const auto neighborIndex = element.mpiIndices[faultSide];
      for (size_t vertexIdx = 0; vertexIdx < Cell::NumVertices; ++vertexIdx) {
        const auto& array3d = ghostMetadataItr->second[neighborIndex].vertices[vertexIdx];
        neighborElemCoords[vertexIdx] = reinterpret_cast<const double (*)[3]>(array3d);
      }
    }

    for (const auto pointId : topology.pointsOf(faceId)) {
      const auto& receiver = outputData_->receivers[topology.representative(pointId)];
      topology.points[pointId].basisFunctions =
          getPlusMinusBasisFunctions(receiver.global.coords, elemCoords, neighborElemCoords);
    }
  }
}

void ReceiverBasedOutputBuilder::initDeviceCollectors(bool elementwise) {
  const auto& faultInfo = meshReader_->getFault();
  const auto& elementsInfo = meshReader_->getElements();

  auto& topology = outputData_->topology;

  std::unordered_map<std::size_t, std::size_t> elementIndices;
  std::unordered_map<std::pair<int, std::size_t>, GhostElement, HashPair<int, std::size_t>>
      elementIndicesGhost;

  // The gather arrays are built per face: a face needs the derivatives of its own element and of
  // its neighbour, no matter how many output points sit on it.
  for (auto& face : topology.faces) {
    const auto elementIndex = faultInfo[face.faultFaceIndex].element;
    const auto& element = elementsInfo[elementIndex];

    if (elementIndices.find(elementIndex) == elementIndices.end()) {
      elementIndices[elementIndex] = elementIndices.size();
    }
    face.deviceDataPlus = elementIndices.at(elementIndex);

    const auto neighborElementIndex = faultInfo[face.faultFaceIndex].neighborElement;
    if (neighborElementIndex >= 0) {
      if (elementIndices.find(neighborElementIndex) == elementIndices.end()) {
        elementIndices[neighborElementIndex] = elementIndices.size();
      }
    } else {
      const auto faultSide = faultInfo[face.faultFaceIndex].side;
      const auto ghostIndex = std::pair<int, std::size_t>(element.neighborRanks[faultSide],
                                                          element.mpiIndices[faultSide]);
      if (elementIndicesGhost.find(ghostIndex) == elementIndicesGhost.end()) {
        const auto index = elementIndicesGhost.size();
        elementIndicesGhost[ghostIndex] =
            GhostElement{std::pair<std::size_t, int>(elementIndex, faultSide), index};
      }
    }
  }

  // the ghost entries are appended behind the local ones, so their offset is only known once all
  // local elements have been seen
  for (auto& face : topology.faces) {
    const auto neighborElementIndex = faultInfo[face.faultFaceIndex].neighborElement;
    if (neighborElementIndex >= 0) {
      face.deviceDataMinus = elementIndices.at(neighborElementIndex);
    } else {
      const auto elementIndex = faultInfo[face.faultFaceIndex].element;
      const auto& element = elementsInfo[elementIndex];
      const auto faultSide = faultInfo[face.faultFaceIndex].side;
      const auto ghostIndex = std::pair<int, std::size_t>(element.neighborRanks[faultSide],
                                                          element.mpiIndices[faultSide]);
      face.deviceDataMinus = elementIndices.size() + elementIndicesGhost.at(ghostIndex).index;
    }
  }

  outputData_->cellCount = elementIndices.size() + elementIndicesGhost.size();

  if constexpr (isDeviceOn()) {
    if (elementwise) {
      // rely on the data that is copied anyways
      // (deviceDataCollector needs to be set to avoid a nullptr call)
      outputData_->deviceDataCollector = std::make_unique<seissol::parallel::DataCollector<real>>(
          std::vector<real*>{}, seissol::kernels::Solver::DerivativesSize, true);
    } else {
      // setup (sparse) index arrays (for receivers)

      std::vector<real*> indexPtrs(outputData_->cellCount);

      for (const auto& [index, arrayIndex] : elementIndices) {
        const auto position = wpBackmap_->get(index);
        indexPtrs[arrayIndex] = wpStorage_->lookup<LTS::DerivativesDevice>(position);
        assert(indexPtrs[arrayIndex] != nullptr);
      }
      for (const auto& [_, ghost] : elementIndicesGhost) {
        const auto neighbor = ghost.data;
        const auto arrayIndex = ghost.index + elementIndices.size();

        const auto position = wpBackmap_->get(neighbor.first);
        indexPtrs[arrayIndex] =
            wpStorage_->lookup<LTS::FaceNeighborsDevice>(position)[neighbor.second];
        assert(indexPtrs[arrayIndex] != nullptr);
      }

      outputData_->deviceDataCollector = std::make_unique<seissol::parallel::DataCollector<real>>(
          indexPtrs, seissol::kernels::Solver::DerivativesSize, useMPIUSM());

      for (const auto& variable : variables_) {
        auto* var = drStorage_->varUntyped(variable, initializer::AllocationPlace::Device);
        const std::size_t elementSize = drStorage_->info(variable).bytes;

        std::vector<void*> dataPointers(topology.faceCount());
        for (std::size_t faceId = 0; faceId < topology.faceCount(); ++faceId) {
          dataPointers[faceId] = reinterpret_cast<uint8_t*>(var) +
                                 elementSize * topology.faces[faceId].position.global;
        }

        const bool hostAccessible = useUSM() && !outputData_->extraRuntime.has_value();
        outputData_->deviceVariables[variable] =
            std::make_unique<seissol::parallel::DataCollectorUntyped>(
                dataPointers, elementSize, hostAccessible);
      }
    }
  }
}

void ReceiverBasedOutputBuilder::initFaultDirections() {
  const auto& faultInfo = meshReader_->getFault();

  for (auto& face : outputData_->topology.faces) {
    auto& faultDirections = face.faultDirections;
    const auto globalIndex = face.faultFaceIndex;

    std::copy_n(&faultInfo[globalIndex].normal[0], 3, faultDirections.faceNormal.begin());
    std::copy_n(&faultInfo[globalIndex].tangent1[0], 3, faultDirections.tangent1.begin());
    std::copy_n(&faultInfo[globalIndex].tangent2[0], 3, faultDirections.tangent2.begin());

    misc::computeStrikeAndDipVectors(faultDirections.faceNormal.data(),
                                     faultDirections.strike.data(),
                                     faultDirections.dip.data());
  }
}

void ReceiverBasedOutputBuilder::initRotationMatrices() {
  using namespace seissol::transformations;
  using RotationMatrixViewT = yateto::DenseTensorView<2, real, unsigned>;

  for (auto& face : outputData_->topology.faces) {
    const auto& faultDirections = face.faultDirections;
    const auto& faceNormal = faultDirections.faceNormal;
    const auto& strike = faultDirections.strike;
    const auto& dip = faultDirections.dip;
    const auto& tangent1 = faultDirections.tangent1;
    const auto& tangent2 = faultDirections.tangent2;

    {
      auto* memorySpace = face.stressGlbToDipStrikeAligned.data();
      RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
      inverseSymmetricTensor2RotationMatrix(
          faceNormal.data(), strike.data(), dip.data(), rotationMatrixView, 0, 0);
    }
    {
      auto* memorySpace = face.stressFaceAlignedToGlb.data();
      RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
      symmetricTensor2RotationMatrix(
          faceNormal.data(), tangent1.data(), tangent2.data(), rotationMatrixView, 0, 0);
    }
    {
      // the face-aligned-to-global direction is not part of the output; it is only needed to
      // obtain its inverse
      std::array<real, seissol::tensor::T::size()> faceAlignedToGlbData{};
      auto faceAlignedToGlb = init::T::view::create(faceAlignedToGlbData.data());
      auto glbToFaceAligned = init::Tinv::view::create(face.glbToFaceAlignedData.data());

      seissol::model::getFaceRotationMatrix(
          faceNormal.data(), tangent1.data(), tangent2.data(), faceAlignedToGlb, glbToFaceAligned);
    }
  }
}

void ReceiverBasedOutputBuilder::initOutputVariables(
    std::array<bool, std::tuple_size_v<DrVarsT>>& outputMask) {
  auto assignMask = [&outputMask](auto& var, int receiverId) {
    var.isActive = outputMask[receiverId];
  };
  misc::forEach(outputData_->vars, assignMask);

  auto allocateVariables = [this](auto& var, int) {
    var.maxCacheLevel = outputData_->maxCacheLevel;
    var.allocateData(this->outputData_->receivers.size());
  };
  misc::forEach(outputData_->vars, allocateVariables);
}

void ReceiverBasedOutputBuilder::initJacobian2dMatrices() {
  const auto& faultInfo = meshReader_->getFault();
  const auto& verticesInfo = meshReader_->getVertices();
  const auto& elementsInfo = meshReader_->getElements();

  for (auto& outputFace : outputData_->topology.faces) {
    const auto& element = elementsInfo[outputFace.elementIndex];
    auto face =
        getGlobalTriangle(static_cast<int>(outputFace.localFaceSideId), element, verticesInfo);

    VrtxCoords xab;
    VrtxCoords xac;
    {
      constexpr size_t X{0};
      constexpr size_t Y{1};
      constexpr size_t Z{2};
      xab[X] = face.point(1)[X] - face.point(0)[X];
      xab[Y] = face.point(1)[Y] - face.point(0)[Y];
      xab[Z] = face.point(1)[Z] - face.point(0)[Z];

      xac[X] = face.point(2)[X] - face.point(0)[X];
      xac[Y] = face.point(2)[Y] - face.point(0)[Y];
      xac[Z] = face.point(2)[Z] - face.point(0)[Z];
    }

    const auto* tangent1 = faultInfo[outputFace.faultFaceIndex].tangent1;
    const auto* tangent2 = faultInfo[outputFace.faultFaceIndex].tangent2;

    Eigen::Matrix<real, 2, 2> matrix;
    matrix(0, 0) = MeshTools::dot(tangent1, xab);
    matrix(0, 1) = MeshTools::dot(tangent2, xab);
    matrix(1, 0) = MeshTools::dot(tangent1, xac);
    matrix(1, 1) = MeshTools::dot(tangent2, xac);
    outputFace.jacobianT2d = matrix.inverse();
  }
}

void ReceiverBasedOutputBuilder::assignNearestInternalGaussianPoints() {
  auto& geoPoints = outputData_->receivers;
  constexpr int NumPoly = ConvergenceOrder - 1;

  for (auto& geoPoint : geoPoints) {
    assert(geoPoint.nearestGpIndex != -1 && "nearestGpIndex must be initialized first");
    if constexpr (Config::DRQuadRule == DRQuadRuleType::Stroud) {
      geoPoint.nearestInternalGpIndex =
          getClosestInternalStroudGp(geoPoint.nearestGpIndex, NumPoly);
    } else {
      geoPoint.nearestInternalGpIndex = geoPoint.nearestGpIndex;
    }
  }
}

void ReceiverBasedOutputBuilder::assignFaultTags() {
  auto& geoPoints = outputData_->receivers;
  const auto& faultInfo = meshReader_->getFault();
  for (auto& geoPoint : geoPoints) {
    geoPoint.faultTag = faultInfo[geoPoint.faultFaceIndex].tag;
  }
}

void ReceiverBasedOutputBuilder::assignFusedIndices() {
  auto& geoPoints = outputData_->receivers;
  for (auto& geoPoint : geoPoints) {
    geoPoint.gpIndex = multisim::NumSimulations * geoPoint.nearestGpIndex + geoPoint.simIndex;
    geoPoint.internalGpIndexFused =
        multisim::NumSimulations * geoPoint.nearestInternalGpIndex + geoPoint.simIndex;
  }
}

} // namespace seissol::dr::output
