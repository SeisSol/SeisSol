// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Constants.h"
#include "DataTypes/ConditionalKey.h"
#include "DataTypes/EncodedConstants.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "Recorders.h"

#include <array>
#include <cstddef>
#include <vector>
#include <yateto.h>

using namespace device;
using namespace seissol::initializer;
using namespace seissol::initializer::recording;

void DynamicRuptureRecorder::record(DynamicRupture::Layer& layer) {
  setUpContext(layer);
  recordDofsTimeEvaluation();
  recordSpaceInterpolation();
}

void DynamicRuptureRecorder::recordDofsTimeEvaluation() {
  real** timeDerivativePlus = currentLayer->var<DynamicRupture::TimeDerivativePlusDevice>();
  real** timeDerivativeMinus = currentLayer->var<DynamicRupture::TimeDerivativeMinusDevice>();
  real* idofsPlus = static_cast<real*>(
      currentLayer->var<DynamicRupture::IdofsPlusOnDevice>(AllocationPlace::Device));
  real* idofsMinus = static_cast<real*>(
      currentLayer->var<DynamicRupture::IdofsMinusOnDevice>(AllocationPlace::Device));

  const auto size = currentLayer->size();
  if (size > 0) {
    std::vector<real*> timeDerivativePlusPtrs(size, nullptr);
    std::vector<real*> timeDerivativeMinusPtrs(size, nullptr);
    std::vector<real*> idofsPlusPtrs(size, nullptr);
    std::vector<real*> idofsMinusPtrs(size, nullptr);

    const size_t idofsSize = tensor::Q::size();
    for (std::size_t faceId = 0; faceId < size; ++faceId) {
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
      currentLayer->var<DynamicRupture::QInterpolatedPlus>(AllocationPlace::Device);
  auto* qInterpolatedMinus =
      currentLayer->var<DynamicRupture::QInterpolatedMinus>(AllocationPlace::Device);

  real* idofsPlus = static_cast<real*>(
      currentLayer->var<DynamicRupture::IdofsPlusOnDevice>(AllocationPlace::Device));
  real* idofsMinus = static_cast<real*>(
      currentLayer->var<DynamicRupture::IdofsMinusOnDevice>(AllocationPlace::Device));

  DRGodunovData* godunovData =
      currentLayer->var<DynamicRupture::GodunovData>(AllocationPlace::Device);
  DRFaceInformation* faceInfo = currentLayer->var<DynamicRupture::FaceInformation>();

  const auto size = currentLayer->size();
  if (size > 0) {
    std::array<std::vector<real*>, *FaceId::Count> qInterpolatedPlusPtr{};
    std::array<std::vector<real*>, *FaceId::Count> idofsPlusPtr{};
    std::array<std::vector<real*>, *FaceId::Count> tInvTPlusPtr{};

    std::array<std::vector<real*>[*FaceId::Count], *FaceId::Count> qInterpolatedMinusPtr {};
    std::array<std::vector<real*>[*FaceId::Count], *FaceId::Count> idofsMinusPtr {};
    std::array<std::vector<real*>[*FaceId::Count], *FaceId::Count> tInvTMinusPtr {};

    const size_t idofsSize = tensor::Q::size();
    for (std::size_t faceId = 0; faceId < size; ++faceId) {
      const auto plusSide = faceInfo[faceId].plusSide;
      qInterpolatedPlusPtr[plusSide].push_back(&qInterpolatedPlus[faceId][0][0]);
      idofsPlusPtr[plusSide].push_back(&idofsPlus[faceId * idofsSize]);
      tInvTPlusPtr[plusSide].push_back((&godunovData[faceId])->dataTinvT);

      const auto minusSide = faceInfo[faceId].minusSide;
      const auto faceRelation = faceInfo[faceId].faceRelation;
      qInterpolatedMinusPtr[minusSide][faceRelation].push_back(&qInterpolatedMinus[faceId][0][0]);
      idofsMinusPtr[minusSide][faceRelation].push_back(&idofsMinus[faceId * idofsSize]);
      tInvTMinusPtr[minusSide][faceRelation].push_back((&godunovData[faceId])->dataTinvT);
    }

    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      if (!qInterpolatedPlusPtr[side].empty()) {
        const ConditionalKey key(*KernelNames::DrSpaceMap, side);
        (*currentDrTable)[key].set(inner_keys::Dr::Id::QInterpolatedPlus,
                                   qInterpolatedPlusPtr[side]);
        (*currentDrTable)[key].set(inner_keys::Dr::Id::IdofsPlus, idofsPlusPtr[side]);
        (*currentDrTable)[key].set(inner_keys::Dr::Id::TinvT, tInvTPlusPtr[side]);
      }
      for (std::size_t faceRelation = 0; faceRelation < Cell::NumFaces; ++faceRelation) {
        if (!qInterpolatedMinusPtr[side][faceRelation].empty()) {
          const ConditionalKey key(*KernelNames::DrSpaceMap, side, faceRelation);
          (*currentDrTable)[key].set(inner_keys::Dr::Id::QInterpolatedMinus,
                                     qInterpolatedMinusPtr[side][faceRelation]);
          (*currentDrTable)[key].set(inner_keys::Dr::Id::IdofsMinus,
                                     idofsMinusPtr[side][faceRelation]);
          (*currentDrTable)[key].set(inner_keys::Dr::Id::TinvT, tInvTMinusPtr[side][faceRelation]);
        }
      }
    }
  }
}
