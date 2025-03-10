// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Recorders.h"
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Tree/Layer.h>
#include <array>
#include <cstddef>
#include <tensor.h>
#include <vector>
#include <yateto.h>

using namespace device;
using namespace seissol::initializer;
using namespace seissol::initializer::recording;

void DynamicRuptureRecorder::record(DynamicRupture& handler, Layer& layer) {
  setUpContext(handler, layer);
  recordDofsTimeEvaluation();
  recordSpaceInterpolation();
}

void DynamicRuptureRecorder::recordDofsTimeEvaluation() {
  real** timeDerivativePlus = currentLayer->var(currentHandler->timeDerivativePlusDevice);
  real** timeDerivativeMinus = currentLayer->var(currentHandler->timeDerivativeMinusDevice);
  real* idofsPlus = static_cast<real*>(currentLayer->getScratchpadMemory(
      currentHandler->idofsPlusOnDevice, AllocationPlace::Device));
  real* idofsMinus = static_cast<real*>(currentLayer->getScratchpadMemory(
      currentHandler->idofsMinusOnDevice, AllocationPlace::Device));

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real*> timeDerivativePlusPtrs(size, nullptr);
    std::vector<real*> timeDerivativeMinusPtrs(size, nullptr);
    std::vector<real*> idofsPlusPtrs(size, nullptr);
    std::vector<real*> idofsMinusPtrs(size, nullptr);

    const size_t idofsSize = tensor::Q::size();
    for (unsigned faceId = 0; faceId < size; ++faceId) {
      timeDerivativePlusPtrs[faceId] = timeDerivativePlus[faceId];
      timeDerivativeMinusPtrs[faceId] = timeDerivativeMinus[faceId];
      idofsPlusPtrs[faceId] = &idofsPlus[faceId * idofsSize];
      idofsMinusPtrs[faceId] = &idofsMinus[faceId * idofsSize];
    }

    const ConditionalKey key(*KernelNames::DrTime);
    checkKey(key);

    (*currentDrTable)[key].set(inner_keys::Dr::Id::DerivativesPlus, timeDerivativePlusPtrs);
    (*currentDrTable)[key].set(inner_keys::Dr::Id::DerivativesMinus, timeDerivativeMinusPtrs);
    (*currentDrTable)[key].set(inner_keys::Dr::Id::IdofsPlus, idofsPlusPtrs);
    (*currentDrTable)[key].set(inner_keys::Dr::Id::IdofsMinus, idofsMinusPtrs);
  }
}

void DynamicRuptureRecorder::recordSpaceInterpolation() {
  auto* qInterpolatedPlus =
      currentLayer->var(currentHandler->qInterpolatedPlus, AllocationPlace::Device);
  auto* qInterpolatedMinus =
      currentLayer->var(currentHandler->qInterpolatedMinus, AllocationPlace::Device);

  real* idofsPlus = static_cast<real*>(currentLayer->getScratchpadMemory(
      currentHandler->idofsPlusOnDevice, AllocationPlace::Device));
  real* idofsMinus = static_cast<real*>(currentLayer->getScratchpadMemory(
      currentHandler->idofsMinusOnDevice, AllocationPlace::Device));

  DRGodunovData* godunovData =
      currentLayer->var(currentHandler->godunovData, AllocationPlace::Device);
  DRFaceInformation* faceInfo = currentLayer->var(currentHandler->faceInformation);

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real*> qInterpolatedPtr{};
    std::vector<real*> idofsPtr{};
    std::vector<real*> tInvTPtr{};
    std::vector<real*> globalPtr{};

    const size_t idofsSize = tensor::Q::size();
    for (unsigned faceId = 0; faceId < size; ++faceId) {
      const auto plusSide = faceInfo[faceId].plusSide;
      qInterpolatedPtr.push_back(&qInterpolatedPlus[faceId][0][0]);
      idofsPtr.push_back(&idofsPlus[faceId * idofsSize]);
      tInvTPtr.push_back((&godunovData[faceId])->TinvT);
      globalPtr.push_back(const_cast<real*>(global->faceToNodalMatrices(plusSide, 0)));

      const auto minusSide = faceInfo[faceId].minusSide;
      const auto faceRelation = faceInfo[faceId].faceRelation;
      qInterpolatedPtr.push_back(&qInterpolatedMinus[faceId][0][0]);
      idofsPtr.push_back(&idofsMinus[faceId * idofsSize]);
      tInvTPtr.push_back((&godunovData[faceId])->TinvT);
      globalPtr.push_back(const_cast<real*>(global->faceToNodalMatrices(minusSide, faceRelation)));
    }

    if (!qInterpolatedPtr.empty()) {
      const ConditionalKey key(*KernelNames::DrSpaceMap);
      (*currentDrTable)[key].set(inner_keys::Dr::Id::QInterpolatedPlus, qInterpolatedPtr);
      (*currentDrTable)[key].set(inner_keys::Dr::Id::IdofsPlus, idofsPtr);
      (*currentDrTable)[key].set(inner_keys::Dr::Id::TinvT, tInvTPtr);
      (*currentDrTable)[key].set(inner_keys::Dr::Id::Global, globalPtr);
    }
  }
}
