// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.h"

#include "Common/Constants.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Output/OutputAux.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "GeneratedCode/tensor.h"
#include "Kernels/Solver.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Helper.h"

#include <memory>
#endif

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
    std::vector<::seissol::initializer::StoragePosition>* faceToLtsMap) {
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

void ReceiverBasedOutputBuilder::initBasisFunctions() {
  const auto& faultInfo = meshReader_->getFault();
  const auto& elementsInfo = meshReader_->getElements();
  const auto& verticesInfo = meshReader_->getVertices();
  const auto& mpiGhostMetadata = meshReader_->getGhostlayerMetadata();

  std::unordered_map<std::size_t, std::size_t> faceIndices;
  std::unordered_map<std::size_t, std::size_t> elementIndices;
  std::unordered_map<std::pair<int, std::size_t>, GhostElement, HashPair<int, std::size_t>>
      elementIndicesGhost;
  std::size_t foundPoints = 0;

  constexpr size_t NumVertices{4};
  for (const auto& point : outputData_->receiverPoints) {
    if (point.isInside) {
      if (faceIndices.find(faceToLtsMap_->at(point.faultFaceIndex).global) == faceIndices.end()) {
        const auto faceIndex = faceIndices.size();
        faceIndices[faceToLtsMap_->at(point.faultFaceIndex).global] = faceIndex;
      }

      ++foundPoints;
      const auto elementIndex = faultInfo[point.faultFaceIndex].element;
      const auto& element = elementsInfo[elementIndex];

      if (elementIndices.find(elementIndex) == elementIndices.end()) {
        const auto index = elementIndices.size();
        elementIndices[elementIndex] = index;
      }

      const auto neighborElementIndex = faultInfo[point.faultFaceIndex].neighborElement;

      const VrtxCoords* elemCoords[NumVertices]{};
      for (size_t vertexIdx = 0; vertexIdx < NumVertices; ++vertexIdx) {
        const auto address = elementsInfo[elementIndex].vertices[vertexIdx];
        elemCoords[vertexIdx] = &(verticesInfo[address].coords);
      }

      const VrtxCoords* neighborElemCoords[NumVertices]{};
      if (neighborElementIndex >= 0) {
        if (elementIndices.find(neighborElementIndex) == elementIndices.end()) {
          const auto index = elementIndices.size();
          elementIndices[neighborElementIndex] = index;
        }
        for (size_t vertexIdx = 0; vertexIdx < NumVertices; ++vertexIdx) {
          const auto address = elementsInfo[neighborElementIndex].vertices[vertexIdx];
          neighborElemCoords[vertexIdx] = &(verticesInfo[address].coords);
        }
      } else {
        const auto faultSide = faultInfo[point.faultFaceIndex].side;
        const auto neighborRank = element.neighborRanks[faultSide];
        const auto& ghostMetadataItr = mpiGhostMetadata.find(neighborRank);
        assert(ghostMetadataItr != mpiGhostMetadata.end());

        const auto neighborIndex = element.mpiIndices[faultSide];

        const auto ghostIndex = std::pair<int, std::size_t>(neighborRank, neighborIndex);
        if (elementIndicesGhost.find(ghostIndex) == elementIndicesGhost.end()) {
          const auto index = elementIndicesGhost.size();
          elementIndicesGhost[ghostIndex] =
              GhostElement{std::pair<std::size_t, int>(elementIndex, faultSide), index};
        }

        for (size_t vertexIdx = 0; vertexIdx < NumVertices; ++vertexIdx) {
          const auto& array3d = ghostMetadataItr->second[neighborIndex].vertices[vertexIdx];
          neighborElemCoords[vertexIdx] = reinterpret_cast<const double (*)[3]>(array3d);
        }
      }

      outputData_->basisFunctions.emplace_back(
          getPlusMinusBasisFunctions(point.global.coords, elemCoords, neighborElemCoords));
    }
  }

  outputData_->cellCount = elementIndices.size() + elementIndicesGhost.size();

#ifdef ACL_DEVICE
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
    indexPtrs[arrayIndex] = wpStorage_->lookup<LTS::FaceNeighborsDevice>(position)[neighbor.second];
    assert(indexPtrs[arrayIndex] != nullptr);
  }

  outputData_->deviceDataCollector = std::make_unique<seissol::parallel::DataCollector<real>>(
      indexPtrs, seissol::kernels::Solver::DerivativesSize, useMPIUSM());

  for (const auto& variable : variables_) {
    auto* var = drStorage_->varUntyped(variable, initializer::AllocationPlace::Device);
    const std::size_t elementSize = drStorage_->info(variable).bytes;

    assert(elementSize % sizeof(real) == 0);

    const std::size_t elementCount = elementSize / sizeof(real);
    std::vector<real*> dataPointers(faceIndices.size());
    for (const auto& [index, arrayIndex] : faceIndices) {
      dataPointers[arrayIndex] = reinterpret_cast<real*>(var) + elementCount * index;
    }

    const bool hostAccessible = useUSM() && !outputData_->extraRuntime.has_value();
    outputData_->deviceVariables[variable] =
        std::make_unique<seissol::parallel::DataCollector<real>>(
            dataPointers, elementCount, hostAccessible);
  }
#endif

  outputData_->deviceDataPlus.resize(foundPoints);
  outputData_->deviceDataMinus.resize(foundPoints);
  outputData_->deviceIndices.resize(foundPoints);
  std::size_t pointCounter = 0;
  for (std::size_t i = 0; i < outputData_->receiverPoints.size(); ++i) {
    const auto& point = outputData_->receiverPoints[i];
    if (point.isInside) {
      const auto elementIndex = faultInfo[point.faultFaceIndex].element;
      const auto& element = elementsInfo[elementIndex];
      outputData_->deviceIndices[pointCounter] =
          faceIndices.at(faceToLtsMap_->at(point.faultFaceIndex).global);

      outputData_->deviceDataPlus[pointCounter] = elementIndices.at(elementIndex);

      const auto neighborElementIndex = faultInfo[point.faultFaceIndex].neighborElement;
      if (neighborElementIndex >= 0) {
        outputData_->deviceDataMinus[pointCounter] = elementIndices.at(neighborElementIndex);
      } else {
        const auto faultSide = faultInfo[point.faultFaceIndex].side;
        const auto neighborRank = element.neighborRanks[faultSide];
        const auto neighborIndex = element.mpiIndices[faultSide];
        outputData_->deviceDataMinus[pointCounter] =
            elementIndices.size() +
            elementIndicesGhost.at(std::pair<int, std::size_t>(neighborRank, neighborIndex)).index;
      }

      ++pointCounter;
    }
  }
}

void ReceiverBasedOutputBuilder::initFaultDirections() {
  const size_t nReceiverPoints = outputData_->receiverPoints.size();
  outputData_->faultDirections.resize(nReceiverPoints);
  const auto& faultInfo = meshReader_->getFault();

  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    const size_t globalIndex = outputData_->receiverPoints[receiverId].faultFaceIndex;

    auto& faceNormal = outputData_->faultDirections[receiverId].faceNormal;
    auto& tangent1 = outputData_->faultDirections[receiverId].tangent1;
    auto& tangent2 = outputData_->faultDirections[receiverId].tangent2;

    std::copy_n(&faultInfo[globalIndex].normal[0], 3, faceNormal.begin());
    std::copy_n(&faultInfo[globalIndex].tangent1[0], 3, tangent1.begin());
    std::copy_n(&faultInfo[globalIndex].tangent2[0], 3, tangent2.begin());

    auto& strike = outputData_->faultDirections[receiverId].strike;
    auto& dip = outputData_->faultDirections[receiverId].dip;
    misc::computeStrikeAndDipVectors(faceNormal.data(), strike.data(), dip.data());
  }
}

void ReceiverBasedOutputBuilder::initRotationMatrices() {
  using namespace seissol::transformations;
  using RotationMatrixViewT = yateto::DenseTensorView<2, real, unsigned>;

  // allocate Rotation Matrices
  // Note: several receiver can share the same rotation matrix
  const size_t nReceiverPoints = outputData_->receiverPoints.size();
  outputData_->stressGlbToDipStrikeAligned.resize(nReceiverPoints);
  outputData_->stressFaceAlignedToGlb.resize(nReceiverPoints);
  outputData_->faceAlignedToGlbData.resize(nReceiverPoints);
  outputData_->glbToFaceAlignedData.resize(nReceiverPoints);

  // init Rotation Matrices
  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    const auto& faceNormal = outputData_->faultDirections[receiverId].faceNormal;
    const auto& strike = outputData_->faultDirections[receiverId].strike;
    const auto& dip = outputData_->faultDirections[receiverId].dip;
    const auto& tangent1 = outputData_->faultDirections[receiverId].tangent1;
    const auto& tangent2 = outputData_->faultDirections[receiverId].tangent2;

    {
      auto* memorySpace = outputData_->stressGlbToDipStrikeAligned[receiverId].data();
      RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
      inverseSymmetricTensor2RotationMatrix(
          faceNormal.data(), strike.data(), dip.data(), rotationMatrixView, 0, 0);
    }
    {
      auto* memorySpace = outputData_->stressFaceAlignedToGlb[receiverId].data();
      RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
      symmetricTensor2RotationMatrix(
          faceNormal.data(), tangent1.data(), tangent2.data(), rotationMatrixView, 0, 0);
    }
    {
      auto faceAlignedToGlb =
          init::T::view::create(outputData_->faceAlignedToGlbData[receiverId].data());
      auto glbToFaceAligned =
          init::Tinv::view::create(outputData_->glbToFaceAlignedData[receiverId].data());

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
    var.allocateData(this->outputData_->receiverPoints.size());
  };
  misc::forEach(outputData_->vars, allocateVariables);
}

void ReceiverBasedOutputBuilder::initJacobian2dMatrices() {
  const auto& faultInfo = meshReader_->getFault();
  const auto& verticesInfo = meshReader_->getVertices();
  const auto& elementsInfo = meshReader_->getElements();

  const size_t nReceiverPoints = outputData_->receiverPoints.size();
  outputData_->jacobianT2d.resize(nReceiverPoints);

  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    const auto side = outputData_->receiverPoints[receiverId].localFaceSideId;
    const auto elementIndex = outputData_->receiverPoints[receiverId].elementIndex;

    assert(elementIndex >= 0);

    const auto& element = elementsInfo[elementIndex];
    auto face = getGlobalTriangle(side, element, verticesInfo);

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

    const auto faultIndex = outputData_->receiverPoints[receiverId].faultFaceIndex;
    const auto* tangent1 = faultInfo[faultIndex].tangent1;
    const auto* tangent2 = faultInfo[faultIndex].tangent2;

    Eigen::Matrix<real, 2, 2> matrix;
    matrix(0, 0) = MeshTools::dot(tangent1, xab);
    matrix(0, 1) = MeshTools::dot(tangent2, xab);
    matrix(1, 0) = MeshTools::dot(tangent1, xac);
    matrix(1, 1) = MeshTools::dot(tangent2, xac);
    outputData_->jacobianT2d[receiverId] = matrix.inverse();
  }
}

void ReceiverBasedOutputBuilder::assignNearestInternalGaussianPoints() {
  auto& geoPoints = outputData_->receiverPoints;
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
  auto& geoPoints = outputData_->receiverPoints;
  const auto& faultInfo = meshReader_->getFault();
  for (auto& geoPoint : geoPoints) {
    geoPoint.faultTag = faultInfo[geoPoint.faultFaceIndex].tag;
  }
}

void ReceiverBasedOutputBuilder::assignFusedIndices() {
  auto& geoPoints = outputData_->receiverPoints;
  for (auto& geoPoint : geoPoints) {
    geoPoint.gpIndex = multisim::NumSimulations * geoPoint.nearestGpIndex + geoPoint.simIndex;
    geoPoint.internalGpIndexFused =
        multisim::NumSimulations * geoPoint.nearestInternalGpIndex + geoPoint.simIndex;
  }
}

} // namespace seissol::dr::output
