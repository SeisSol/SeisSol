#include "Recorders.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;


void DynamicRuptureRecorder::record(DynamicRupture &handler, Layer &layer) {
  setUpContext(handler, layer);
  recordDofsTimeEvaluation();
  recordSpaceInterpolation();
}


void DynamicRuptureRecorder::recordDofsTimeEvaluation() {
  real** timeDerivativePlus = currentLayer->var(currentHandler->timeDerivativePlus);
  real** timeDerivativeMinus = currentLayer->var(currentHandler->timeDerivativeMinus);
  real* idofsPlus = static_cast<real *>(currentLayer->getScratchpadMemory(currentHandler->idofsPlusOnDevice));
  real* idofsMinus = static_cast<real *>(currentLayer->getScratchpadMemory(currentHandler->idofsMinusOnDevice));

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::vector<real *> timeDerivativePlusPtrs(size, nullptr);
    std::vector<real *> timeDerivativeMinusPtrs(size, nullptr);
    std::vector<real *> idofsPlusPtrs(size, nullptr);
    std::vector<real *> idofsMinusPtrs(size, nullptr);

    const size_t idofsSize = tensor::Q::size();
    for (unsigned faceId = 0; faceId < size; ++faceId) {
      timeDerivativePlusPtrs[faceId] = timeDerivativePlus[faceId];
      timeDerivativeMinusPtrs[faceId] = timeDerivativeMinus[faceId];
      idofsPlusPtrs[faceId] = &idofsPlus[faceId * idofsSize];
      idofsMinusPtrs[faceId] = &idofsMinus[faceId * idofsSize];
    }

    ConditionalKey key(*KernelNames::DrTime);
    checkKey(key);

    (*currentTable)[key].content[*EntityId::DrDerivativesPlus] = new BatchPointers(timeDerivativePlusPtrs);
    (*currentTable)[key].content[*EntityId::DrDerivativesMinus] = new BatchPointers(timeDerivativeMinusPtrs);
    (*currentTable)[key].content[*EntityId::DrIdofsPlus] = new BatchPointers(idofsPlusPtrs);
    (*currentTable)[key].content[*EntityId::DrIdofsMinus] = new BatchPointers(idofsMinusPtrs);
  }
}


void DynamicRuptureRecorder::recordSpaceInterpolation() {
  auto* qInterpolatedPlus = currentLayer->var(currentHandler->qInterpolatedPlus);
  auto* qInterpolatedMinus = currentLayer->var(currentHandler->qInterpolatedMinus);

  real* idofsPlus = static_cast<real *>(currentLayer->getScratchpadMemory(currentHandler->idofsPlusOnDevice));
  real* idofsMinus = static_cast<real *>(currentLayer->getScratchpadMemory(currentHandler->idofsMinusOnDevice));

  DRGodunovData* godunovData = currentLayer->var(currentHandler->godunovData);
  DRFaceInformation* faceInfo = currentLayer->var(currentHandler->faceInformation);

  const auto size = currentLayer->getNumberOfCells();
  if (size > 0) {
    std::array<std::vector<real *>, *FaceId::Count> qInterpolatedPlusPtr {};
    std::array<std::vector<real *>, *FaceId::Count> idofsPlusPtr {};
    std::array<std::vector<real *>, *FaceId::Count> TinvTPlusPtr {};

    std::array<std::vector<real *>[*FaceId::Count], *FaceId::Count> qInterpolatedMinusPtr {};
    std::array<std::vector<real *>[*FaceId::Count], *FaceId::Count> idofsMinusPtr {};
    std::array<std::vector<real *>[*FaceId::Count], *FaceId::Count> TinvTMinusPtr {};

    const size_t idofsSize = tensor::Q::size();
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
        (*currentTable)[key].content[*EntityId::DrQInterpolatedPlus] = new BatchPointers(qInterpolatedPlusPtr[side]);
        (*currentTable)[key].content[*EntityId::DrIdofsPlus] = new BatchPointers(idofsPlusPtr[side]);
        (*currentTable)[key].content[*EntityId::DrTinvT] = new BatchPointers(TinvTPlusPtr[side]);
      }
      for (unsigned faceRelation = 0; faceRelation < 4; ++faceRelation) {
        if (!qInterpolatedMinusPtr[side][faceRelation].empty()) {
          ConditionalKey key(*KernelNames::DrSpaceMap, side, faceRelation);
          (*currentTable)[key].content[*EntityId::DrQInterpolatedMinus] = new BatchPointers(qInterpolatedMinusPtr[side][faceRelation]);
          (*currentTable)[key].content[*EntityId::DrIdofsMinus] = new BatchPointers(idofsMinusPtr[side][faceRelation]);
          (*currentTable)[key].content[*EntityId::DrTinvT] = new BatchPointers(TinvTMinusPtr[side][faceRelation]);
        }
      }
    }
  }
}
