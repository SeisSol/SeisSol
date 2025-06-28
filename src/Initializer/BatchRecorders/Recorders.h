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
#include "Memory/Tree/Layer.h"
#include "utils/logger.h"
#include <vector>

namespace seissol::initializer::recording {
class AbstractRecorder {
  public:
  virtual ~AbstractRecorder() = default;

  virtual void record(Layer& layer) = 0;

  protected:
  void checkKey(const ConditionalKey& key) {
    if (currentTable->find(key) != currentTable->end()) {
      logError()
          << "Table key conflict detected. Problems with hashing in batch recording subsystem";
    }
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
  Layer* currentLayer{nullptr};
};

class CompositeRecorder : public AbstractRecorder {
  public:
  ~CompositeRecorder() override {
    for (auto* recorder : concreteRecorders) {
      delete recorder;
    }
  }

  void record(Layer& layer) override {
    for (auto* recorder : concreteRecorders) {
      recorder->record(layer);
    }
  }

  void addRecorder(AbstractRecorder* recorder) { concreteRecorders.push_back(recorder); }

  void removeRecorder(size_t recorderIndex) {
    if (recorderIndex < concreteRecorders.size()) {
      concreteRecorders.erase(concreteRecorders.begin() + recorderIndex);
    }
  }

  private:
  std::vector<AbstractRecorder*> concreteRecorders;
};

class LocalIntegrationRecorder : public AbstractRecorder {
  public:
  void record(Layer& layer) override;

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

class NeighIntegrationRecorder : public AbstractRecorder {
  public:
  void record(Layer& layer) override;

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

class PlasticityRecorder : public AbstractRecorder {
  public:
  void setUpContext(Layer& layer) { AbstractRecorder::setUpContext(layer); }

  void record(Layer& layer) override;
};

class DynamicRuptureRecorder : public AbstractRecorder {
  public:
  void record(Layer& layer) override;

  private:
  void setUpContext(Layer& layer) { AbstractRecorder::setUpContext(layer); }
  void recordDofsTimeEvaluation();
  void recordSpaceInterpolation();
  std::unordered_map<real*, real*> idofsAddressRegistry;
};

} // namespace seissol::initializer::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
