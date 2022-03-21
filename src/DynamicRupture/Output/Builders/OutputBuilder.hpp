#ifndef SEISSOL_DR_OUTPUT_BUILDER_HPP
#define SEISSOL_DR_OUTPUT_BUILDER_HPP

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.hpp"
#include "DynamicRupture/Output/OutputAux.hpp"
#include "Geometry/MeshReader.h"
#include "Initializer/InputAux.hpp"
#include "Model/common.hpp"
#include "Numerical_aux/Transformation.h"
#include "Parallel/MPI.h"

namespace seissol::dr::output {
class OutputBuilder {
  public:
  virtual ~OutputBuilder() {
    auto deallocateVars = [](auto& var, int) { var.releaseData(); };
    misc::forEach(outputData.vars, deallocateVars);
  }

  virtual void init() = 0;
  virtual void initTimeCaching() = 0;

  void setMeshReader(const MeshReader* reader) {
    meshReader = reader;
    localRank = MPI::mpi.rank();
  }

  void initBasisFunctions() {
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    for (const auto& point : outputData.receiverPoints) {
      if (point.isInside) {
        auto elementIndex = faultInfo[point.faultFaceIndex].element;
        auto neighborElementIndex = faultInfo[point.faultFaceIndex].neighborElement;

        const VrtxCoords* elemCoords[4]{};
        const VrtxCoords* neighborElemCoords[4]{};

        for (int i = 0; i < 4; ++i) {
          elemCoords[i] = &(verticesInfo[elementsInfo[elementIndex].vertices[i]].coords);
          neighborElemCoords[i] =
              &(verticesInfo[elementsInfo[neighborElementIndex].vertices[i]].coords);
        }

        outputData.basisFunctions.emplace_back(
            getPlusMinusBasisFunctions(point.global.coords, elemCoords, neighborElemCoords));
      }
    }
  }

  void initFaultDirections() {
    size_t nReceiverPoints = outputData.receiverPoints.size();
    outputData.faultDirections.resize(nReceiverPoints);
    const auto& faultInfo = meshReader->getFault();

    for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
      size_t globalIndex = outputData.receiverPoints[receiverId].faultFaceIndex;

      outputData.faultDirections[receiverId].faceNormal = faultInfo[globalIndex].normal;
      outputData.faultDirections[receiverId].tangent1 = faultInfo[globalIndex].tangent1;
      outputData.faultDirections[receiverId].tangent2 = faultInfo[globalIndex].tangent2;
      misc::computeStrikeAndDipVectors(outputData.faultDirections[receiverId].faceNormal,
                                       outputData.faultDirections[receiverId].strike,
                                       outputData.faultDirections[receiverId].dip);
    }
  }

  void initRotationMatrices() {
    using namespace seissol::transformations;
    using RotationMatrixViewT = yateto::DenseTensorView<2, real, unsigned>;

    // allocate Rotation Matrices
    // Note: several receiver can share the same rotation matrix
    size_t nReceiverPoints = outputData.receiverPoints.size();
    outputData.stressGlbToDipStrikeAligned.resize(nReceiverPoints);
    outputData.stressFaceAlignedToGlb.resize(nReceiverPoints);
    outputData.faceAlignedToGlbData.resize(nReceiverPoints);
    outputData.glbToFaceAlignedData.resize(nReceiverPoints);

    // init Rotation Matrices
    for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
      const auto* const faceNormal = outputData.faultDirections[receiverId].faceNormal;
      auto* const strike = outputData.faultDirections[receiverId].strike;
      auto* const dip = outputData.faultDirections[receiverId].dip;
      const auto* const tangent1 = outputData.faultDirections[receiverId].tangent1;
      const auto* const tangent2 = outputData.faultDirections[receiverId].tangent2;

      {
        auto* memorySpace = outputData.stressGlbToDipStrikeAligned[receiverId].data();
        RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
        inverseSymmetricTensor2RotationMatrix(faceNormal, strike, dip, rotationMatrixView, 0, 0);
      }
      {
        auto* memorySpace = outputData.stressFaceAlignedToGlb[receiverId].data();
        RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
        symmetricTensor2RotationMatrix(faceNormal, tangent1, tangent2, rotationMatrixView, 0, 0);
      }
      {
        auto faceAlignedToGlb =
            init::T::view::create(outputData.faceAlignedToGlbData[receiverId].data());
        auto glbToFaceAligned =
            init::Tinv::view::create(outputData.glbToFaceAlignedData[receiverId].data());

        seissol::model::getFaceRotationMatrix(
            faceNormal, tangent1, tangent2, faceAlignedToGlb, glbToFaceAligned);
      }
    }
  }

  void initOutputVariables(std::array<bool, std::tuple_size<DrVarsT>::value>& outputMask) {
    auto assignMask = [&outputMask](auto& var, int receiverId) {
      var.isActive = outputMask[receiverId];
    };
    misc::forEach(outputData.vars, assignMask);

    auto allocateVariables = [this](auto& var, int) {
      var.maxCacheLevel = outputData.maxCacheLevel;
      var.allocateData(this->outputData.receiverPoints.size());
    };
    misc::forEach(outputData.vars, allocateVariables);
  }

  void initJacobian2dMatrices() {
    const auto& faultInfo = meshReader->getFault();
    const auto& verticesInfo = meshReader->getVertices();
    const auto& elementsInfo = meshReader->getElements();

    size_t nReceiverPoints = outputData.receiverPoints.size();
    outputData.jacobianT2d.resize(nReceiverPoints);

    for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
      assert(outputData.receiverPoints[receiverId].elementIndex != -1);

      auto side = outputData.receiverPoints[receiverId].localFaceSideId;
      auto elementIndex = outputData.receiverPoints[receiverId].elementIndex;

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

      auto faultIndex = outputData.receiverPoints[receiverId].faultFaceIndex;
      auto* tangent1 = faultInfo[faultIndex].tangent1;
      auto* tangent2 = faultInfo[faultIndex].tangent2;

      Eigen::Matrix<real, 2, 2> matrix;
      matrix(0, 0) = MeshTools::dot(tangent1, xab);
      matrix(0, 1) = MeshTools::dot(tangent2, xab);
      matrix(1, 0) = MeshTools::dot(tangent1, xac);
      matrix(1, 1) = MeshTools::dot(tangent2, xac);
      outputData.jacobianT2d[receiverId] = matrix.inverse();
    }
  }

  protected:
  const MeshReader* meshReader{};
  OutputData outputData;
  int localRank{-1};
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_BUILDER_HPP
