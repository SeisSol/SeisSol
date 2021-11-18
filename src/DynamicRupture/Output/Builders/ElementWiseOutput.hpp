#ifndef SEISSOL_DRELEMENTWISEOUTPUT_HPP
#define SEISSOL_DRELEMENTWISEOUTPUT_HPP

#include <DynamicRupture/Math.h>
#include "DynamicRupture/Output/DataTypes.hpp"
#include "DynamicRupture/Output/FaultRefiner/FaultRefiner.hpp"
#include "Geometry/MeshReader.h"
#include "DynamicRupture/Output/OutputAux.hpp"
#include "Parallel/MPI.h"

namespace seissol::dr::output {
class Base;

struct FaultGeomParamsT {
  int numSubTriangles{};
  int numSubElements{};
  size_t numSides{};
};

class ElementWiseOutput {
  public:
  friend Base;

  ~ElementWiseOutput() {
    auto deallocateVars = [](auto& var, int) { var.releaseData(); };
    aux::forEach(outputData.vars, deallocateVars);
  }

  void setParams(const ElementwiseFaultParamsT& params, const MeshReader* reader) {
    elementwiseParams = params;
    meshReader = reader;
    localRank = MPI::mpi.rank();
  }

  void init() {
    initReceiverLocations();
    assignNearestGaussianPoints(outputData.receiverPoints);
    initTimeCaching();
    initOutputVariables();
    initFaultDirections();
    initRotationMatrices();
    initBasisFunctions();
    outputData.isActive = true;
  }

  void initTimeCaching() {
    outputData.maxCacheLevel = ElementWiseOutput::maxAllowedCacheLevel;
    outputData.currentCacheLevel = 0;
  }

  void initOutputVariables() {
    auto assignMask = [this](auto& var, int index) {
      var.isActive = this->elementwiseParams.outputMask[index];
    };
    aux::forEach(outputData.vars, assignMask);

    auto allocateVariables = [this](auto& var, int) {
      var.maxCacheLevel = outputData.maxCacheLevel;
      var.allocateData(this->outputData.receiverPoints.size());
    };
    aux::forEach(outputData.vars, allocateVariables);
  }

  void initReceiverLocations() {
    std::unique_ptr<FaultRefinerInterface> faultRefiner{nullptr};
    faultRefiner = getRefiner(elementwiseParams.refinementStrategy);

    geomParam.numSides = meshReader->getFault().size();
    geomParam.numSubTriangles = faultRefiner->getNumSubTriangles();
    geomParam.numSubElements = std::pow(geomParam.numSubTriangles, elementwiseParams.refinement);

    logInfo(localRank) << "CPP: Initialising Fault output. Refinement strategy: "
                       << elementwiseParams.refinementStrategy
                       << " Number of sub-triangles: " << geomParam.numSubTriangles;

    // get arrays of elements and vertices from the mesher
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    // iterate through each fault side
    for (size_t faceIndex = 0; faceIndex < geomParam.numSides; ++faceIndex) {

      // get a Global Element ID for the current fault face
      auto elementIndex = faultInfo[faceIndex].element;
      if (elementIndex > 0) {

        // store coords of vertices of the current ELEMENT
        std::array<const double*, 4> elementVerticesCoords{};
        for (int ElementVertexId = 0; ElementVertexId < 4; ++ElementVertexId) {
          auto globalVertexId = elementsInfo[elementIndex].vertices[ElementVertexId];
          elementVerticesCoords[ElementVertexId] = verticesInfo[globalVertexId].coords;
        }

        auto localFaceSideId = faultInfo[faceIndex].side;

        // init reference coordinates of the fault face
        ExtTriangle referenceFace = getReferenceFace(localFaceSideId);

        // init global coordinates of the fault face
        ExtTriangle globalFace = getReferenceFace(localFaceSideId);
        for (int faceVertexId = 0; faceVertexId < 3; ++faceVertexId) {
          auto elementVertexId = getElementVertexId(localFaceSideId, faceVertexId);
          auto globalVertexId = elementsInfo[elementIndex].vertices[elementVertexId];

          for (int coordId = 0; coordId < 3; ++coordId) {
            globalFace.points[faceVertexId][coordId] = verticesInfo[globalVertexId].coords[coordId];
          }
        }

        faultRefiner->refineAndAccumulate(
            elementwiseParams.refinement, faceIndex, localFaceSideId, referenceFace, globalFace);
      }
    }

    // retrieve all receivers from a fault face refiner
    outputData.receiverPoints = faultRefiner->moveAllReceiverPoints();
    faultRefiner.reset(nullptr);

    nOutPoints = outputData.receiverPoints.size();
  }

  void initFaultDirections() {
    outputData.faultDirections.resize(geomParam.numSides);
    const auto& faultInfo = meshReader->getFault();

    for (size_t index = 0; index < geomParam.numSides; ++index) {
      outputData.faultDirections[index].faceNormal = faultInfo[index].normal;
      outputData.faultDirections[index].tangent1 = faultInfo[index].tangent1;
      outputData.faultDirections[index].tangent2 = faultInfo[index].tangent2;
      computeStrikeAndDipVectors(outputData.faultDirections[index].faceNormal,
                                 outputData.faultDirections[index].strike,
                                 outputData.faultDirections[index].dip);
    }
  }

  void initRotationMatrices() {
    using namespace seissol::transformations;
    using RotationMatrixViewT = yateto::DenseTensorView<2, real, unsigned>;

    // allocate Rotation Matrices
    // Note: several receiver can share the same rotation matrix
    outputData.rotationMatrices.resize(geomParam.numSides);

    // init Rotation Matrices
    for (size_t index = 0; index < geomParam.numSides; ++index) {
      const auto faceNormal = outputData.faultDirections[index].faceNormal;
      const auto strike = outputData.faultDirections[index].strike;
      const auto dip = outputData.faultDirections[index].dip;

      computeStrikeAndDipVectors(faceNormal, strike, dip);

      std::vector<real> rotationMatrix(36, 0.0);
      RotationMatrixViewT rotationMatrixView(const_cast<real*>(rotationMatrix.data()), {6, 6});

      symmetricTensor2RotationMatrix(faceNormal, strike, dip, rotationMatrixView, 0, 0);
      outputData.rotationMatrices[index] = std::move(rotationMatrix);
    }
  }

  void initBasisFunctions() {
    const auto& faultInfo = meshReader->getFault();
    const auto& elementsInfo = meshReader->getElements();
    const auto& verticesInfo = meshReader->getVertices();

    for (const auto& point : outputData.receiverPoints) {
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

  void initConstrains() {
    /*
    for (const auto& Point: m_ReceiverPoints) {
      const auto& RotationMatrix = m_RotationMatrices[Point.FaultFaceIndex];

    }
    */
  }

  void evaluateInitialStressInFaultCS() {
    // Compute initialStressInFaultCS
  }

  inline const static size_t maxAllowedCacheLevel = 1;

  private:
  ElementwiseFaultParamsT elementwiseParams;
  FaultGeomParamsT geomParam;

  const MeshReader* meshReader{nullptr};
  size_t nOutPoints{};
  int localRank{-1};
  OutputData outputData{};
};
} // namespace seissol::dr::output
#endif // SEISSOL_DRELEMENTWISEOUTPUT_HPP
