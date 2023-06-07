#include "DynamicRupture/Output/Builders/ReceiverBasedOutputBuilder.hpp"
#include "Model/common.hpp"

namespace seissol::dr::output {
void ReceiverBasedOutputBuilder::setMeshReader(const seissol::geometry::MeshReader* reader) {
  meshReader = reader;
  localRank = MPI::mpi.rank();
}

void ReceiverBasedOutputBuilder::initBasisFunctions() {
  const auto& faultInfo = meshReader->getFault();
  const auto& elementsInfo = meshReader->getElements();
  const auto& verticesInfo = meshReader->getVertices();
  const auto& mpiGhostMetadata = meshReader->getGhostlayerMetadata();

  constexpr size_t numVertices{4};
  for (const auto& point : outputData->receiverPoints) {
    if (point.isInside) {
      const auto elementIndex = faultInfo[point.faultFaceIndex].element;
      const auto& element = elementsInfo[elementIndex];

      const auto neighborElementIndex = faultInfo[point.faultFaceIndex].neighborElement;

      const VrtxCoords* elemCoords[numVertices]{};
      for (size_t vertexIdx = 0; vertexIdx < numVertices; ++vertexIdx) {
        const auto address = elementsInfo[elementIndex].vertices[vertexIdx];
        elemCoords[vertexIdx] = &(verticesInfo[address].coords);
      }

      const VrtxCoords* neighborElemCoords[numVertices]{};
      if (neighborElementIndex >= 0) {
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

      seissol::model::getFaceRotationMatrix<seissol::model::Material_t>(
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
  constexpr int numPoly = ConvergenceOrder - 1;

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
