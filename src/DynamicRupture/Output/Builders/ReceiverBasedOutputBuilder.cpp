// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.h"
#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Output/OutputAux.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Lut.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <init.h>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Parallel/DataCollector.h"
#include "Parallel/Helper.h"
#include <Memory/Tree/Layer.h>
#include <memory>
#include <tensor.h>
#endif

namespace seissol::dr::output {
void ReceiverBasedOutputBuilder::setMeshReader(const seissol::geometry::MeshReader* reader) {
  meshReader = reader;
  localRank = MPI::mpi.rank();
}

void ReceiverBasedOutputBuilder::setLtsData(seissol::initializer::LTSTree* userWpTree,
                                            seissol::initializer::LTS* userWpDescr,
                                            seissol::initializer::Lut* userWpLut,
                                            seissol::initializer::LTSTree* userDrTree,
                                            seissol::initializer::DynamicRupture* userDrDescr) {
  wpTree = userWpTree;
  wpDescr = userWpDescr;
  wpLut = userWpLut;
  drTree = userDrTree;
  drDescr = userDrDescr;
}

void ReceiverBasedOutputBuilder::setVariableList(const std::vector<std::size_t>& variables) {
  this->variables = variables;
}

void ReceiverBasedOutputBuilder::setFaceToLtsMap(std::vector<std::size_t>* faceToLtsMap) {
  this->faceToLtsMap = faceToLtsMap;
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
    seed ^= hasher2(data.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};
} // namespace

void ReceiverBasedOutputBuilder::initBasisFunctions() {
  const auto& faultInfo = meshReader->getFault();
  const auto& elementsInfo = meshReader->getElements();
  const auto& verticesInfo = meshReader->getVertices();
  const auto& mpiGhostMetadata = meshReader->getGhostlayerMetadata();

  std::unordered_map<std::size_t, std::size_t> faceIndices;
  std::unordered_map<std::size_t, std::size_t> elementIndices;
  std::unordered_map<std::pair<int, std::size_t>, GhostElement, HashPair<int, std::size_t>>
      elementIndicesGhost;
  std::size_t foundPoints = 0;

  constexpr size_t NumVertices{4};
  for (const auto& point : outputData->receiverPoints) {
    if (point.isInside) {
      if (faceIndices.find(faceToLtsMap->at(point.faultFaceIndex)) == faceIndices.end()) {
        const auto faceIndex = faceIndices.size();
        faceIndices[faceToLtsMap->at(point.faultFaceIndex)] = faceIndex;
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
          auto* data = const_cast<double*>(array3d);
          neighborElemCoords[vertexIdx] = reinterpret_cast<double(*)[3]>(data);
        }
      }

      outputData->basisFunctions.emplace_back(
          getPlusMinusBasisFunctions(point.global.coords, elemCoords, neighborElemCoords));
    }
  }

  outputData->cellCount = elementIndices.size() + elementIndicesGhost.size();

#ifdef ACL_DEVICE
  std::vector<real*> indexPtrs(outputData->cellCount);

  for (const auto& [index, arrayIndex] : elementIndices) {
    indexPtrs[arrayIndex] = wpLut->lookup(wpDescr->derivativesDevice, index);
    assert(indexPtrs[arrayIndex] != nullptr);
  }
  for (const auto& [_, ghost] : elementIndicesGhost) {
    const auto neighbor = ghost.data;
    const auto arrayIndex = ghost.index + elementIndices.size();
    indexPtrs[arrayIndex] =
        wpLut->lookup(wpDescr->faceNeighborsDevice, neighbor.first)[neighbor.second];
    assert(indexPtrs[arrayIndex] != nullptr);
  }

  outputData->deviceDataCollector = std::make_unique<seissol::parallel::DataCollector>(
      indexPtrs, seissol::tensor::Q::size(), useMPIUSM());

  for (const auto& variable : variables) {
    auto* var = drTree->varUntyped(variable, initializer::AllocationPlace::Device);
    const std::size_t elementSize = drTree->info(variable).elemsize;

    assert(elementSize % sizeof(real) == 0);

    const std::size_t elementCount = elementSize / sizeof(real);
    std::vector<real*> dataPointers(faceIndices.size());
    for (const auto& [index, arrayIndex] : faceIndices) {
      dataPointers[arrayIndex] = reinterpret_cast<real*>(var) + elementCount * index;
    }
    outputData->deviceVariables[variable] =
        std::make_unique<seissol::parallel::DataCollector>(dataPointers, elementCount, useUSM());
  }
#endif

  outputData->deviceDataPlus.resize(foundPoints);
  outputData->deviceDataMinus.resize(foundPoints);
  outputData->deviceIndices.resize(foundPoints);
  std::size_t pointCounter = 0;
  for (std::size_t i = 0; i < outputData->receiverPoints.size(); ++i) {
    const auto& point = outputData->receiverPoints[i];
    if (point.isInside) {
      const auto elementIndex = faultInfo[point.faultFaceIndex].element;
      const auto& element = elementsInfo[elementIndex];
      outputData->deviceIndices[pointCounter] =
          faceIndices.at(faceToLtsMap->at(point.faultFaceIndex));

      outputData->deviceDataPlus[pointCounter] = elementIndices.at(elementIndex);

      const auto neighborElementIndex = faultInfo[point.faultFaceIndex].neighborElement;
      if (neighborElementIndex >= 0) {
        outputData->deviceDataMinus[pointCounter] = elementIndices.at(neighborElementIndex);
      } else {
        const auto faultSide = faultInfo[point.faultFaceIndex].side;
        const auto neighborRank = element.neighborRanks[faultSide];
        const auto neighborIndex = element.mpiIndices[faultSide];
        outputData->deviceDataMinus[pointCounter] =
            elementIndices.size() +
            elementIndicesGhost.at(std::pair<int, std::size_t>(neighborRank, neighborIndex)).index;
      }

      ++pointCounter;
    }
  }
}

void ReceiverBasedOutputBuilder::initFaultDirections() {
  const size_t nReceiverPoints = outputData->receiverPoints.size();
  outputData->faultDirections.resize(nReceiverPoints);
  const auto& faultInfo = meshReader->getFault();

  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    const size_t globalIndex = outputData->receiverPoints[receiverId].faultFaceIndex;

    auto& faceNormal = outputData->faultDirections[receiverId].faceNormal;
    auto& tangent1 = outputData->faultDirections[receiverId].tangent1;
    auto& tangent2 = outputData->faultDirections[receiverId].tangent2;

    std::copy_n(&faultInfo[globalIndex].normal[0], 3, faceNormal.begin());
    std::copy_n(&faultInfo[globalIndex].tangent1[0], 3, tangent1.begin());
    std::copy_n(&faultInfo[globalIndex].tangent2[0], 3, tangent2.begin());

    auto& strike = outputData->faultDirections[receiverId].strike;
    auto& dip = outputData->faultDirections[receiverId].dip;
    misc::computeStrikeAndDipVectors(faceNormal.data(), strike.data(), dip.data());
  }
}

void ReceiverBasedOutputBuilder::initRotationMatrices() {
  using namespace seissol::transformations;
  using RotationMatrixViewT = yateto::DenseTensorView<2, real, unsigned>;

  // allocate Rotation Matrices
  // Note: several receiver can share the same rotation matrix
  const size_t nReceiverPoints = outputData->receiverPoints.size();
  outputData->stressGlbToDipStrikeAligned.resize(nReceiverPoints);
  outputData->stressFaceAlignedToGlb.resize(nReceiverPoints);
  outputData->faceAlignedToGlbData.resize(nReceiverPoints);
  outputData->glbToFaceAlignedData.resize(nReceiverPoints);

  // init Rotation Matrices
  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    const auto& faceNormal = outputData->faultDirections[receiverId].faceNormal;
    const auto& strike = outputData->faultDirections[receiverId].strike;
    const auto& dip = outputData->faultDirections[receiverId].dip;
    const auto& tangent1 = outputData->faultDirections[receiverId].tangent1;
    const auto& tangent2 = outputData->faultDirections[receiverId].tangent2;

    {
      auto* memorySpace = outputData->stressGlbToDipStrikeAligned[receiverId].data();
      RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
      inverseSymmetricTensor2RotationMatrix(
          faceNormal.data(), strike.data(), dip.data(), rotationMatrixView, 0, 0);
    }
    {
      auto* memorySpace = outputData->stressFaceAlignedToGlb[receiverId].data();
      RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
      symmetricTensor2RotationMatrix(
          faceNormal.data(), tangent1.data(), tangent2.data(), rotationMatrixView, 0, 0);
    }
    {
      auto faceAlignedToGlb =
          init::T::view::create(outputData->faceAlignedToGlbData[receiverId].data());
      auto glbToFaceAligned =
          init::Tinv::view::create(outputData->glbToFaceAlignedData[receiverId].data());

      seissol::model::getFaceRotationMatrix(
          faceNormal.data(), tangent1.data(), tangent2.data(), faceAlignedToGlb, glbToFaceAligned);
    }
  }
}

void ReceiverBasedOutputBuilder::initOutputVariables(
    std::array<bool, std::tuple_size<DrVarsT>::value>& outputMask) {
  auto assignMask = [&outputMask](auto& var, int receiverId) {
    var.isActive = outputMask[receiverId];
  };
  misc::forEach(outputData->vars, assignMask);

  auto allocateVariables = [this](auto& var, int) {
    var.maxCacheLevel = outputData->maxCacheLevel;
    var.allocateData(this->outputData->receiverPoints.size());
  };
  misc::forEach(outputData->vars, allocateVariables);
}

void ReceiverBasedOutputBuilder::initJacobian2dMatrices() {
  const auto& faultInfo = meshReader->getFault();
  const auto& verticesInfo = meshReader->getVertices();
  const auto& elementsInfo = meshReader->getElements();

  const size_t nReceiverPoints = outputData->receiverPoints.size();
  outputData->jacobianT2d.resize(nReceiverPoints);

  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    const auto side = outputData->receiverPoints[receiverId].localFaceSideId;
    const auto elementIndex = outputData->receiverPoints[receiverId].elementIndex;

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

    const auto faultIndex = outputData->receiverPoints[receiverId].faultFaceIndex;
    const auto* tangent1 = faultInfo[faultIndex].tangent1;
    const auto* tangent2 = faultInfo[faultIndex].tangent2;

    Eigen::Matrix<real, 2, 2> matrix;
    matrix(0, 0) = MeshTools::dot(tangent1, xab);
    matrix(0, 1) = MeshTools::dot(tangent2, xab);
    matrix(1, 0) = MeshTools::dot(tangent1, xac);
    matrix(1, 1) = MeshTools::dot(tangent2, xac);
    outputData->jacobianT2d[receiverId] = matrix.inverse();
  }
}

void ReceiverBasedOutputBuilder::assignNearestInternalGaussianPoints() {
  auto& geoPoints = outputData->receiverPoints;
  constexpr int NumPoly = ConvergenceOrder - 1;

  for (auto& geoPoint : geoPoints) {
    assert(geoPoint.nearestGpIndex != -1 && "nearestGpIndex must be initialized first");
#ifdef stroud
    geoPoint.nearestInternalGpIndex = getClosestInternalStroudGp(geoPoint.nearestGpIndex, NumPoly);
#else
    geoPoint.nearestInternalGpIndex = geoPoint.nearestGpIndex;
#endif
  }
}

void ReceiverBasedOutputBuilder::assignFaultTags() {
  auto& geoPoints = outputData->receiverPoints;
  const auto& faultInfo = meshReader->getFault();
  for (auto& geoPoint : geoPoints) {
    geoPoint.faultTag = faultInfo[geoPoint.faultFaceIndex].tag;
  }
}

} // namespace seissol::dr::output
