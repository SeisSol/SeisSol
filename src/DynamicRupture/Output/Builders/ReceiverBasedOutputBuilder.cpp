#include "DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.hpp"

namespace seissol::dr::output {
void ReceiverBasedOutputBuilder::setMeshReader(const MeshReader* reader) {
  meshReader = reader;
  localRank = MPI::mpi.rank();
}

void ReceiverBasedOutputBuilder::initBasisFunctions() {
  const auto& faultInfo = meshReader->getFault();
  const auto& elementsInfo = meshReader->getElements();
  const auto& verticesInfo = meshReader->getVertices();
  const auto& mpiNeighborVertices = meshReader->getMPINeighborVertices();

  constexpr size_t numVertices{4};
  for (const auto& point : outputData->receiverPoints) {
    if (point.isInside) {
      auto elementIndex = faultInfo[point.faultFaceIndex].element;
      auto& element = elementsInfo[elementIndex];

      auto neighborElementIndex = faultInfo[point.faultFaceIndex].neighborElement;

      const VrtxCoords* elemCoords[numVertices]{};
      for (size_t vertexIdx = 0; vertexIdx < numVertices; ++vertexIdx) {
        auto address = elementsInfo[elementIndex].vertices[vertexIdx];
        elemCoords[vertexIdx] = &(verticesInfo[address].coords);
      }

      const VrtxCoords* neighborElemCoords[numVertices]{};
      if (neighborElementIndex >= 0) {
        for (size_t vertexIdx = 0; vertexIdx < numVertices; ++vertexIdx) {
          auto address = elementsInfo[neighborElementIndex].vertices[vertexIdx];
          neighborElemCoords[vertexIdx] = &(verticesInfo[address].coords);
        }
      } else {
        auto faultSide = faultInfo[point.faultFaceIndex].side;
        auto neighborRank = element.neighborRanks[faultSide];
        const auto& neighborVerticesItr = mpiNeighborVertices.find(neighborRank);
        assert(neighborVerticesItr != mpiNeighborVertices.end());

        auto neighborIndex = element.mpiIndices[faultSide];
        for (size_t vertexIdx = 0; vertexIdx < numVertices; ++vertexIdx) {
          const auto& array3d = neighborVerticesItr->second[neighborIndex][vertexIdx];
          auto* data = const_cast<double*>(array3d.data());
          neighborElemCoords[vertexIdx] = reinterpret_cast<double(*)[3]>(data);
        }
      }

      outputData->basisFunctions.emplace_back(
          getPlusMinusBasisFunctions(point.global.coords, elemCoords, neighborElemCoords));
    }
  }
}

void ReceiverBasedOutputBuilder::initFaultDirections() {
  size_t nReceiverPoints = outputData->receiverPoints.size();
  outputData->faultDirections.resize(nReceiverPoints);
  const auto& faultInfo = meshReader->getFault();

  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    size_t globalIndex = outputData->receiverPoints[receiverId].faultFaceIndex;

    outputData->faultDirections[receiverId].faceNormal = faultInfo[globalIndex].normal;
    outputData->faultDirections[receiverId].tangent1 = faultInfo[globalIndex].tangent1;
    outputData->faultDirections[receiverId].tangent2 = faultInfo[globalIndex].tangent2;
    misc::computeStrikeAndDipVectors(outputData->faultDirections[receiverId].faceNormal,
                                     outputData->faultDirections[receiverId].strike,
                                     outputData->faultDirections[receiverId].dip);
  }
}

void ReceiverBasedOutputBuilder::initRotationMatrices() {
  using namespace seissol::transformations;
  using RotationMatrixViewT = yateto::DenseTensorView<2, real, unsigned>;

  // allocate Rotation Matrices
  // Note: several receiver can share the same rotation matrix
  size_t nReceiverPoints = outputData->receiverPoints.size();
  outputData->stressGlbToDipStrikeAligned.resize(nReceiverPoints);
  outputData->stressFaceAlignedToGlb.resize(nReceiverPoints);
  outputData->faceAlignedToGlbData.resize(nReceiverPoints);
  outputData->glbToFaceAlignedData.resize(nReceiverPoints);

  // init Rotation Matrices
  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    const auto* const faceNormal = outputData->faultDirections[receiverId].faceNormal;
    const auto* const strike = outputData->faultDirections[receiverId].strike;
    const auto* const dip = outputData->faultDirections[receiverId].dip;
    const auto* const tangent1 = outputData->faultDirections[receiverId].tangent1;
    const auto* const tangent2 = outputData->faultDirections[receiverId].tangent2;

    {
      auto* memorySpace = outputData->stressGlbToDipStrikeAligned[receiverId].data();
      RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
      inverseSymmetricTensor2RotationMatrix(faceNormal, strike, dip, rotationMatrixView, 0, 0);
    }
    {
      auto* memorySpace = outputData->stressFaceAlignedToGlb[receiverId].data();
      RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
      symmetricTensor2RotationMatrix(faceNormal, tangent1, tangent2, rotationMatrixView, 0, 0);
    }
    {
      auto faceAlignedToGlb =
          init::T::view::create(outputData->faceAlignedToGlbData[receiverId].data());
      auto glbToFaceAligned =
          init::Tinv::view::create(outputData->glbToFaceAlignedData[receiverId].data());

      seissol::model::getFaceRotationMatrix(
          faceNormal, tangent1, tangent2, faceAlignedToGlb, glbToFaceAligned);
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

  size_t nReceiverPoints = outputData->receiverPoints.size();
  outputData->jacobianT2d.resize(nReceiverPoints);

  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    auto side = outputData->receiverPoints[receiverId].localFaceSideId;
    auto elementIndex = outputData->receiverPoints[receiverId].elementIndex;

    assert(elementIndex >= 0);

    const auto& element = elementsInfo[elementIndex];
    auto face = getGlobalTriangle(side, element, verticesInfo);

    VrtxCoords xab, xac;
    {
      constexpr size_t x{0}, y{1}, z{2};
      xab[x] = face.p2[x] - face.p1[x];
      xab[y] = face.p2[y] - face.p1[y];
      xab[z] = face.p2[z] - face.p1[z];

      xac[x] = face.p3[x] - face.p1[x];
      xac[y] = face.p3[y] - face.p1[y];
      xac[z] = face.p3[z] - face.p1[z];
    }

    auto faultIndex = outputData->receiverPoints[receiverId].faultFaceIndex;
    auto* tangent1 = faultInfo[faultIndex].tangent1;
    auto* tangent2 = faultInfo[faultIndex].tangent2;

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

#ifndef stroud
  logWarning() << "internal gaussian points are equal to "
               << "the nearest gaussian points of the fault receivers";
#endif

  for (auto& geoPoint : geoPoints) {
    assert(geoPoint.nearestGpIndex != -1 && "nearestGpIndex must be initialized first");
#ifdef stroud
    geoPoint.nearestInternalGpIndex = getClosestInternalStroudGp(geoPoint.nearestGpIndex, numPoly);
#else
    geoPoint.nearestInternalGpIndex = geoPoint.nearestGpIndex;
#endif
  }
}
} // namespace seissol::dr::output
