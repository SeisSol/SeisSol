#include "Recorders.h"
#include <Kernels/Interface.hpp>
#include <Common/cellconfigconv.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

template<typename Config>
void DynamicRuptureRecorder<Config>::record(DynamicRupture<Config>& handler, Layer& layer) {
  setUpContext(handler, layer);
  recordDofsTimeEvaluation();
  recordSpaceInterpolation();
}

template<typename Config>
void DynamicRuptureRecorder<Config>::recordDofsTimeEvaluation() {
  RealT** timeDerivativePlus = currentLayer->var(currentHandler->timeDerivativePlus);
  RealT** timeDerivativeMinus = currentLayer->var(currentHandler->timeDerivativeMinus);
  RealT* idofsPlus =
      static_cast<RealT*>(currentLayer->getScratchpadMemory(currentHandler->idofsPlusOnDevice));
  RealT* idofsMinus =
      static_cast<RealT*>(currentLayer->getScratchpadMemory(currentHandler->idofsMinusOnDevice));

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<RealT*> timeDerivativePlusPtrs(size, nullptr);
    std::vector<RealT*> timeDerivativeMinusPtrs(size, nullptr);
    std::vector<RealT*> idofsPlusPtrs(size, nullptr);
    std::vector<RealT*> idofsMinusPtrs(size, nullptr);

    const size_t idofsSize = Yateto<Config>::Tensor::Q::size();
    for (unsigned faceId = 0; faceId < size; ++faceId) {
      timeDerivativePlusPtrs[faceId] = timeDerivativePlus[faceId];
      timeDerivativeMinusPtrs[faceId] = timeDerivativeMinus[faceId];
      idofsPlusPtrs[faceId] = &idofsPlus[faceId * idofsSize];
      idofsMinusPtrs[faceId] = &idofsMinus[faceId * idofsSize];
    }

    ConditionalKey key(*KernelNames::DrTime);
    checkKey(key);

    (*currentDrTable)[key].set(inner_keys::Dr::Id::DerivativesPlus, timeDerivativePlusPtrs);
    (*currentDrTable)[key].set(inner_keys::Dr::Id::DerivativesMinus, timeDerivativeMinusPtrs);
    (*currentDrTable)[key].set(inner_keys::Dr::Id::IdofsPlus, idofsPlusPtrs);
    (*currentDrTable)[key].set(inner_keys::Dr::Id::IdofsMinus, idofsMinusPtrs);
  }
}

template<typename Config>
void DynamicRuptureRecorder<Config>::recordSpaceInterpolation() {
  auto* qInterpolatedPlus = currentLayer->var(currentHandler->qInterpolatedPlus);
  auto* qInterpolatedMinus = currentLayer->var(currentHandler->qInterpolatedMinus);

  RealT* idofsPlus =
      static_cast<RealT*>(currentLayer->getScratchpadMemory(currentHandler->idofsPlusOnDevice));
  RealT* idofsMinus =
      static_cast<RealT*>(currentLayer->getScratchpadMemory(currentHandler->idofsMinusOnDevice));

  DRGodunovData* godunovData = currentLayer->var(currentHandler->godunovData);
  DRFaceInformation* faceInfo = currentLayer->var(currentHandler->faceInformation);

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::array<std::vector<RealT*>, *FaceId::Count> qInterpolatedPlusPtr{};
    std::array<std::vector<RealT*>, *FaceId::Count> idofsPlusPtr{};
    std::array<std::vector<RealT*>, *FaceId::Count> TinvTPlusPtr{};

    std::array<std::vector<RealT*>[*FaceId::Count], *FaceId::Count> qInterpolatedMinusPtr {};
    std::array<std::vector<RealT*>[*FaceId::Count], *FaceId::Count> idofsMinusPtr {};
    std::array<std::vector<RealT*>[*FaceId::Count], *FaceId::Count> TinvTMinusPtr {};

    const size_t idofsSize = Yateto<Config>::Tensor::Q::size();
    for (unsigned faceId = 0; faceId < size; ++faceId) {
      const auto plusSide = faceInfo[faceId].plusSide;
      qInterpolatedPlusPtr[plusSide].push_back(&qInterpolatedPlus[faceId][0][0]);
      idofsPlusPtr[plusSide].push_back(&idofsPlus[faceId * idofsSize]);
      TinvTPlusPtr[plusSide].push_back((&godunovData[faceId])->TinvT);

      const auto minusSide = faceInfo[faceId].minusSide;
      const auto faceRelation = faceInfo[faceId].faceRelation;
      qInterpolatedMinusPtr[minusSide][faceRelation].push_back(&qInterpolatedMinus[faceId][0][0]);
      idofsMinusPtr[minusSide][faceRelation].push_back(&idofsMinus[faceId * idofsSize]);
      TinvTMinusPtr[minusSide][faceRelation].push_back((&godunovData[faceId])->TinvT);
    }

    for (unsigned side = 0; side < 4; ++side) {
      if (!qInterpolatedPlusPtr[side].empty()) {
        ConditionalKey key(*KernelNames::DrSpaceMap, side);
        (*currentDrTable)[key].set(inner_keys::Dr::Id::QInterpolatedPlus,
                                   qInterpolatedPlusPtr[side]);
        (*currentDrTable)[key].set(inner_keys::Dr::Id::IdofsPlus, idofsPlusPtr[side]);
        (*currentDrTable)[key].set(inner_keys::Dr::Id::TinvT, TinvTPlusPtr[side]);
      }
      for (unsigned faceRelation = 0; faceRelation < 4; ++faceRelation) {
        if (!qInterpolatedMinusPtr[side][faceRelation].empty()) {
          ConditionalKey key(*KernelNames::DrSpaceMap, side, faceRelation);
          (*currentDrTable)[key].set(inner_keys::Dr::Id::QInterpolatedMinus,
                                     qInterpolatedMinusPtr[side][faceRelation]);
          (*currentDrTable)[key].set(inner_keys::Dr::Id::IdofsMinus,
                                     idofsMinusPtr[side][faceRelation]);
          (*currentDrTable)[key].set(inner_keys::Dr::Id::TinvT, TinvTMinusPtr[side][faceRelation]);
        }
      }
    }
  }
}

const DeclareForAllConfigs<DynamicRuptureRecorder> declDynamicRuptureRecorder;
