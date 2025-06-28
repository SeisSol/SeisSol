// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_

#include "DataTypes/ConditionalTable.h"
#include "Kernels/Interface.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "utils/logger.h"
#include <vector>

namespace seissol::initializer::recording {
template <typename LtsT>
class AbstractRecorder {
  public:
  virtual ~AbstractRecorder() = default;

  virtual void record(const LtsT& lts, Layer& layer) = 0;

  protected:
  void checkKey(const ConditionalKey& key) {
    if (currentTable->find(key) != currentTable->end()) {
      logError()
          << "Table key conflict detected. Problems with hashing in batch recording subsystem";
    }
  }

  void setUpContext(const LtsT& handler, Layer& layer) {
    currentTable = &(layer.getConditionalTable<inner_keys::Wp>());
    currentDrTable = &(layer.getConditionalTable<inner_keys::Dr>());
    currentMaterialTable = &(layer.getConditionalTable<inner_keys::Material>());
    currentIndicesTable = &(layer.getConditionalTable<inner_keys::Indices>());
    currentHandler = &(handler);
    currentLayer = &(layer);
  }

  void setUpContext(Layer& layer) {
    currentTable = &(layer.getConditionalTable<inner_keys::Wp>());
    currentDrTable = &(layer.getConditionalTable<inner_keys::Dr>());
    currentMaterialTable = &(layer.getConditionalTable<inner_keys::Material>());
    currentIndicesTable = &(layer.getConditionalTable<inner_keys::Indices>());
    currentLayer = &(layer);
  }

  ConditionalPointersToRealsTable* currentTable{nullptr};
  DrConditionalPointersToRealsTable* currentDrTable{nullptr};
  ConditionalMaterialTable* currentMaterialTable{nullptr};
  ConditionalIndicesTable* currentIndicesTable{nullptr};
  const LtsT* currentHandler{nullptr};
  Layer* currentLayer{nullptr};
};

template <typename LtsT>
class CompositeRecorder : public AbstractRecorder<LtsT> {
  public:
  ~CompositeRecorder() override {
    for (auto recorder : concreteRecorders) {
      delete recorder;
    }
  }

  void record(const LtsT& handler, Layer& layer) override {
    for (auto recorder : concreteRecorders) {
      recorder->record(handler, layer);
    }
  }

  void addRecorder(AbstractRecorder<LtsT>* recorder) { concreteRecorders.push_back(recorder); }

  void removeRecorder(size_t recorderIndex) {
    if (recorderIndex < concreteRecorders.size()) {
      concreteRecorders.erase(concreteRecorders.begin() + recorderIndex);
    }
  }

  private:
  std::vector<AbstractRecorder<LtsT>*> concreteRecorders{};
};

class LocalIntegrationRecorder : public AbstractRecorder<seissol::LTS> {
  public:
  void record(const seissol::LTS& lts, Layer& layer) override;

  private:
  void setUpContext(Layer& layer) {
    integratedDofsAddressCounter = 0;
    derivativesAddressCounter = 0;
    AbstractRecorder::setUpContext(layer);
  }

  void recordTimeAndVolumeIntegrals();
  void recordFreeSurfaceGravityBc();
  void recordDirichletBc();
  void recordAnalyticalBc(Layer& layer);
  void recordLocalFluxIntegral();
  void recordDisplacements();

  std::unordered_map<size_t, real*> idofsAddressRegistry;
  std::vector<real*> dQPtrs;

  size_t integratedDofsAddressCounter{0};
  size_t derivativesAddressCounter{0};
};

class NeighIntegrationRecorder : public AbstractRecorder<seissol::LTS> {
  public:
  void record(const seissol::LTS& lts, Layer& layer) override;

  private:
  void setUpContext(Layer& layer) {
    integratedDofsAddressCounter = 0;
    AbstractRecorder::setUpContext(layer);
  }
  void recordDofsTimeEvaluation();
  void recordNeighborFluxIntegrals();
  std::unordered_map<real*, real*> idofsAddressRegistry;
  size_t integratedDofsAddressCounter{0};
};

class PlasticityRecorder : public AbstractRecorder<seissol::LTS> {
  public:
  void setUpContext(Layer& layer) { AbstractRecorder::setUpContext(layer); }

  void record(const seissol::LTS& lts, Layer& layer) override;
};

class DynamicRuptureRecorder : public AbstractRecorder<seissol::DynamicRupture> {
  public:
  void record(const DynamicRupture& handler, Layer& layer) override;

  private:
  void setUpContext(const DynamicRupture& handler, Layer& layer) {
    AbstractRecorder::setUpContext(handler, layer);
  }
  void recordDofsTimeEvaluation();
  void recordSpaceInterpolation();
  std::unordered_map<real*, real*> idofsAddressRegistry;
};

} // namespace seissol::initializer::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
