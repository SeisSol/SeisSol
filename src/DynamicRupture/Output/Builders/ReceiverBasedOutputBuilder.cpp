#include "DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.hpp"
#include <unordered_map>

namespace seissol::dr::output {
void ReceiverBasedOutputBuilder::setMeshReader(const seissol::geometry::MeshReader* reader) {
  meshReader = reader;
  localRank = MPI::mpi.rank();
}

void ReceiverBasedOutputBuilder::setLtsData(seissol::initializer::LTSTree* userWpTree,
                                            seissol::initializer::LTS* userWpDescr,
                                            seissol::initializer::Lut* userWpLut) {
  wpTree = userWpTree;
  wpDescr = userWpDescr;
  wpLut = userWpLut;
}

namespace {
struct GhostElement {
  std::pair<std::size_t, int> data;
  std::size_t index;
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

  std::unordered_map<std::size_t, std::size_t> elementIndices;
  std::unordered_map<std::pair<int, std::size_t>, GhostElement, HashPair<int, std::size_t>>
      elementIndicesGhost;
  std::size_t foundPoints = 0;

  constexpr size_t numVertices{4};
  for (const auto& point : outputData->receiverPoints) {
    if (point.isInside) {
      ++foundPoints;
      const auto elementIndex = faultInfo[point.faultFaceIndex].element;
      const auto& element = elementsInfo[elementIndex];

      if (elementIndices.find(elementIndex) == elementIndices.end()) {
        const auto index = elementIndices.size();
        elementIndices[elementIndex] = index;
      }

      const auto neighborElementIndex = faultInfo[point.faultFaceIndex].neighborElement;

      const VrtxCoords* elemCoords[numVertices]{};
      for (size_t vertexIdx = 0; vertexIdx < numVertices; ++vertexIdx) {
        const auto address = elementsInfo[elementIndex].vertices[vertexIdx];
        elemCoords[vertexIdx] = &(verticesInfo[address].coords);
      }

      const VrtxCoords* neighborElemCoords[numVertices]{};
      if (neighborElementIndex >= 0) {
        if (elementIndices.find(neighborElementIndex) == elementIndices.end()) {
          const auto index = elementIndices.size();
          elementIndices[neighborElementIndex] = index;
        }
        for (size_t vertexIdx = 0; vertexIdx < numVertices; ++vertexIdx) {
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

        for (size_t vertexIdx = 0; vertexIdx < numVertices; ++vertexIdx) {
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
  real** preDofPtr =
      reinterpret_cast<real**>(device::DeviceInstance::getInstance().api->allocPinnedMem(
          sizeof(real*) * outputData->cellCount));
  outputData->deviceDataPtr =
      reinterpret_cast<real**>(device::DeviceInstance::getInstance().api->allocGlobMem(
          sizeof(real*) * outputData->cellCount));

  for (const auto& [index, arrayIndex] : elementIndices) {
    preDofPtr[arrayIndex] = wpLut->lookup(wpDescr->derivatives, index);
    assert(preDofPtr[arrayIndex] != nullptr);
  }
  for (const auto& [_, ghost] : elementIndicesGhost) {
    const auto neighbor = ghost.data;
    const auto arrayIndex = ghost.index + elementIndices.size();
    preDofPtr[arrayIndex] = wpLut->lookup(wpDescr->faceNeighbors, neighbor.first)[neighbor.second];
    assert(preDofPtr[arrayIndex] != nullptr);
  }

  device::DeviceInstance::getInstance().api->copyTo(
      outputData->deviceDataPtr, preDofPtr, sizeof(real*) * outputData->cellCount);
  device::DeviceInstance::getInstance().api->freePinnedMem(preDofPtr);
#endif

  outputData->deviceDataPlus.resize(foundPoints);
  outputData->deviceDataMinus.resize(foundPoints);
  std::size_t pointCounter = 0;
  for (std::size_t i = 0; i < outputData->receiverPoints.size(); ++i) {
    const auto& point = outputData->receiverPoints[i];
    if (point.isInside) {
      const auto elementIndex = faultInfo[point.faultFaceIndex].element;
      const auto& element = elementsInfo[elementIndex];
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

    VrtxCoords xab, xac;
    {
      constexpr size_t x{0}, y{1}, z{2};
      xab[x] = face.point(1)[x] - face.point(0)[x];
      xab[y] = face.point(1)[y] - face.point(0)[y];
      xab[z] = face.point(1)[z] - face.point(0)[z];

      xac[x] = face.point(2)[x] - face.point(0)[x];
      xac[y] = face.point(2)[y] - face.point(0)[y];
      xac[z] = face.point(2)[z] - face.point(0)[z];
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
  constexpr int numPoly = CONVERGENCE_ORDER - 1;

  for (auto& geoPoint : geoPoints) {
    assert(geoPoint.nearestGpIndex != -1 && "nearestGpIndex must be initialized first");
#ifdef stroud
    geoPoint.nearestInternalGpIndex = getClosestInternalStroudGp(geoPoint.nearestGpIndex, numPoly);
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
