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
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include "Solver/MultipleSimulations.h"

#include <Common/ConfigHelper.h>
#include <Common/Typedefs.h>
#include <Eigen/Core>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Descriptor/LTS.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "GeneratedCode/tensor.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Helper.h"

#include <memory>
#endif

namespace seissol::dr::output {
void ReceiverBasedOutputBuilder::setMeshReader(const seissol::geometry::MeshReader* reader) {
  meshReader = reader;
  localRank = Mpi::mpi.rank();
}

void ReceiverBasedOutputBuilder::setLtsData(LTS::Storage& userWpStorage,
                                            LTS::Backmap& userWpBackmap,
                                            DynamicRupture::Storage& userDrStorage) {
  wpStorage = &userWpStorage;
  wpBackmap = &userWpBackmap;
  drStorage = &userDrStorage;
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
    seed ^= hasher2(data.second) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
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

  if (outputData->transformData.empty()) {
    outputData->transformData.resize(outputData->receiverPoints.size());
  }

  std::unordered_set<std::size_t> configIds;

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

      const std::size_t configId = element.configId;
      configIds.insert(configId);

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
          neighborElemCoords[vertexIdx] = reinterpret_cast<double (*)[3]>(data);
        }
      }

      std::visit(
          [&](auto cfg) {
            using Cfg = decltype(cfg);
            if (outputData->transformData[foundPoints - 1].index() != configId) {
              outputData->transformData[foundPoints - 1].emplace<TransformData<Cfg>>();
            }

            std::get<TransformData<Cfg>>(outputData->transformData[foundPoints - 1])
                .basisFunctions = getPlusMinusBasisFunctions<Real<Cfg>>(
                point.global.coords, elemCoords, neighborElemCoords, cfg);
          },
          ConfigVariantList[configId]);
    }
  }

  outputData->cellCount = elementIndices.size() + elementIndicesGhost.size();

  if (configIds.size() != 1) {
    logError() << "Error while setting up receiver data.";
  }

#ifdef ACL_DEVICE
  std::visit(
      [&](auto cfg) {
        using Cfg = decltype(cfg);
        using real = Real<Cfg>;
        std::vector<real*> indexPtrs(outputData->cellCount);

        for (const auto& [index, arrayIndex] : elementIndices) {
          const auto position = wpBackmap->get(index);
          indexPtrs[arrayIndex] = wpStorage->lookup<LTS::DerivativesDevice>(cfg, position);
          assert(indexPtrs[arrayIndex] != nullptr);
        }
        for (const auto& [_, ghost] : elementIndicesGhost) {
          const auto neighbor = ghost.data;
          const auto arrayIndex = ghost.index + elementIndices.size();

          const auto position = wpBackmap->get(neighbor.first);
          indexPtrs[arrayIndex] = reinterpret_cast<real*>(
              wpStorage->lookup<LTS::FaceNeighborsDevice>(position)[neighbor.second]);
          assert(indexPtrs[arrayIndex] != nullptr);
        }

        outputData->deviceDataCollector =
            std::make_unique<DataCollectorVariant>(seissol::parallel::DataCollector<real>(
                indexPtrs, seissol::tensor::Q<Cfg>::size(), useMPIUSM()));

        for (const auto& variable : variables) {
          auto* var = drStorage->varUntyped(variable, initializer::AllocationPlace::Device);
          const std::size_t elementSize = drStorage->info(variable).bytes;

          assert(elementSize % sizeof(real) == 0);

          const std::size_t elementCount = elementSize / sizeof(real);
          std::vector<real*> dataPointers(faceIndices.size());
          for (const auto& [index, arrayIndex] : faceIndices) {
            dataPointers[arrayIndex] = reinterpret_cast<real*>(var) + elementCount * index;
          }

          const bool hostAccessible = useUSM() && !outputData->extraRuntime.has_value();
          outputData->deviceVariables[variable] = std::make_unique<DataCollectorVariant>(
              seissol::parallel::DataCollector<real>(dataPointers, elementCount, hostAccessible));
        }
      },
      ConfigVariantList[*configIds.begin()]);
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

  // allocate Rotation Matrices
  // Note: several receiver can share the same rotation matrix
  const size_t nReceiverPoints = outputData->receiverPoints.size();
  if (outputData->transformData.empty()) {
    outputData->transformData.resize(nReceiverPoints);
  }

  // init Rotation Matrices
  for (size_t receiverId = 0; receiverId < nReceiverPoints; ++receiverId) {
    const size_t configId =
        meshReader->getElements()[outputData->receiverPoints[receiverId].elementIndex].configId;

    std::visit(
        [&](auto cfg) {
          using Cfg = decltype(cfg);
          using real = Real<Cfg>;

          if (outputData->transformData[receiverId].index() != configId) {
            outputData->transformData[receiverId].emplace<TransformData<Cfg>>();
          }
          auto& transformData = std::get<TransformData<Cfg>>(outputData->transformData[receiverId]);

          using RotationMatrixViewT = yateto::DenseTensorView<2, real, unsigned>;
          const auto& faceNormal = outputData->faultDirections[receiverId].faceNormal;
          const auto& strike = outputData->faultDirections[receiverId].strike;
          const auto& dip = outputData->faultDirections[receiverId].dip;
          const auto& tangent1 = outputData->faultDirections[receiverId].tangent1;
          const auto& tangent2 = outputData->faultDirections[receiverId].tangent2;

          {
            auto* memorySpace = transformData.stressGlbToDipStrikeAligned.data();
            RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
            inverseSymmetricTensor2RotationMatrix(
                faceNormal.data(), strike.data(), dip.data(), rotationMatrixView, 0, 0);
          }
          {
            auto* memorySpace = transformData.stressFaceAlignedToGlb.data();
            RotationMatrixViewT rotationMatrixView(memorySpace, {6, 6});
            symmetricTensor2RotationMatrix(
                faceNormal.data(), tangent1.data(), tangent2.data(), rotationMatrixView, 0, 0);
          }
          {
            auto faceAlignedToGlb =
                init::T<Cfg>::view::create(transformData.faceAlignedToGlbData.data());
            auto glbToFaceAligned =
                init::Tinv<Cfg>::view::create(transformData.glbToFaceAlignedData.data());

            seissol::model::getFaceRotationMatrix<Cfg>(faceNormal.data(),
                                                       tangent1.data(),
                                                       tangent2.data(),
                                                       faceAlignedToGlb,
                                                       glbToFaceAligned);
          }
        },
        ConfigVariantList[configId]);
  }
}

void ReceiverBasedOutputBuilder::initOutputVariables(
    std::array<bool, std::tuple_size_v<DrVarsT>>& outputMask) {
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
  if (outputData->transformData.empty()) {
    outputData->transformData.resize(nReceiverPoints);
  }

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

    const size_t configId =
        meshReader->getElements()[outputData->receiverPoints[receiverId].elementIndex].configId;
    std::visit(
        [&](auto cfg) {
          using Cfg = decltype(cfg);
          using real = Real<Cfg>;

          if (outputData->transformData[receiverId].index() != configId) {
            outputData->transformData[receiverId].emplace<TransformData<Cfg>>();
          }

          Eigen::Matrix<real, 2, 2> matrix;
          matrix(0, 0) = MeshTools::dot(tangent1, xab);
          matrix(0, 1) = MeshTools::dot(tangent2, xab);
          matrix(1, 0) = MeshTools::dot(tangent1, xac);
          matrix(1, 1) = MeshTools::dot(tangent2, xac);
          std::get<TransformData<Cfg>>(outputData->transformData[receiverId]).jacobianT2d =
              matrix.inverse();
        },
        ConfigVariantList[configId]);
  }
}

void ReceiverBasedOutputBuilder::assignNearestInternalGaussianPoints() {
  auto& geoPoints = outputData->receiverPoints;

  for (auto& geoPoint : geoPoints) {
    const std::size_t configId = meshReader->getElements()[geoPoint.elementIndex].configId;
    std::visit(
        [&](auto cfg) {
          using Cfg = decltype(cfg);
          constexpr int NumPoly = Cfg::ConvergenceOrder - 1;
          assert(geoPoint.nearestGpIndex != -1 && "nearestGpIndex must be initialized first");
          if constexpr (Cfg::DRQuadRule == DRQuadRuleType::Stroud) {
            geoPoint.nearestInternalGpIndex =
                getClosestInternalStroudGp(geoPoint.nearestGpIndex, NumPoly);
          } else {
            geoPoint.nearestInternalGpIndex = geoPoint.nearestGpIndex;
          }
        },
        ConfigVariantList[configId]);
  }
}

void ReceiverBasedOutputBuilder::assignFaultTags() {
  auto& geoPoints = outputData->receiverPoints;
  const auto& faultInfo = meshReader->getFault();
  for (auto& geoPoint : geoPoints) {
    geoPoint.faultTag = faultInfo[geoPoint.faultFaceIndex].tag;
  }
}

void ReceiverBasedOutputBuilder::assignFusedIndices() {
  auto& geoPoints = outputData->receiverPoints;
  for (auto& geoPoint : geoPoints) {
    const auto& element = meshReader->getElements().at(geoPoint.elementIndex);
    std::size_t simcount = 1;
    std::visit([&](auto cfg) { simcount = decltype(cfg)::NumSimulations; },
               ConfigVariantList[element.configId]);

    geoPoint.gpIndex = simcount * geoPoint.nearestGpIndex + geoPoint.simIndex;
    geoPoint.internalGpIndexFused = simcount * geoPoint.nearestInternalGpIndex + geoPoint.simIndex;
  }
}

} // namespace seissol::dr::output
