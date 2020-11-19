#ifndef SEISSOL_RECORDERS_H
#define SEISSOL_RECORDERS_H

#include "DataTypes/ConditionalTable.hpp"
#include "utils/logger.h"
#include <Initializer/LTS.h>
#include <Initializer/tree/Layer.hpp>
#include <Kernels/Interface.hpp>
#include <vector>

namespace seissol {
namespace initializers {
namespace recording {


class AbstractRecorder {
public:
  virtual ~AbstractRecorder() = default;

  virtual void record(LTS &handler, Layer &layer) = 0;

protected:
  void checkKey(const ConditionalKey &key) {
    if (currentTable->find(key) != currentTable->end()) {
      logError()
          << "Table key conflict detected. Problems with hashing in batch recording subsystem";
    }
  }

  void setUpContext(LTS &handler, Layer &layer) {
    currentTable = &(layer.getCondBatchTable());
    currentHandler = &(handler);
    currentLayer = &(layer);

    idofsAddressCounter = 0;
    derivativesAddressCounter = 0;
  }

  ConditionalBatchTableT *currentTable{nullptr};
  LTS *currentHandler{nullptr};
  Layer *currentLayer{nullptr};

  size_t idofsAddressCounter{0};
  size_t derivativesAddressCounter{0};
};


class CompositeRecorder : public AbstractRecorder {
public:
  ~CompositeRecorder() override {
    for (auto recorder : concreteRecorders)
      delete recorder;
  }

  void record(LTS &handler, Layer &layer) override {
    for (auto recorder : concreteRecorders) {
      recorder->record(handler, layer);
    }
  }

  void addRecorder(AbstractRecorder *recorder) {
    concreteRecorders.push_back(recorder);
  }

  void removeRecorder(size_t recorderIndex) {
    if (recorderIndex < concreteRecorders.size()) {
      concreteRecorders.erase(concreteRecorders.begin() + recorderIndex);
    }
  }

private:
  std::vector<AbstractRecorder *> concreteRecorders{};
};


class LocalIntegrationRecorder : public AbstractRecorder {
public:
  void record(LTS &handler, Layer &layer) override;

private:
  void setUpContext(LTS &handler, Layer &layer, kernels::LocalData::Loader &loader) {
    currentLoader = &loader;
    AbstractRecorder::setUpContext(handler, layer);
  }

  kernels::LocalData::Loader *currentLoader{nullptr};
  void recordTimeAndVolumeIntegrals();
  void recordLocalFluxIntegral();
  void recordDisplacements();
  std::unordered_map<size_t, real *> idofsAddressRegistry{};
};


class NeighIntegrationRecorder : public AbstractRecorder {
public:
  void record(LTS &handler, Layer &layer) override;
private:
  void setUpContext(LTS &handler, Layer &layer, kernels::NeighborData::Loader &loader) {
    currentLoader = &loader;
    AbstractRecorder::setUpContext(handler, layer);
  }
  void recordDofsTimeEvaluation();
  void recordNeighbourFluxIntegrals();
  kernels::NeighborData::Loader *currentLoader{nullptr};
  std::unordered_map<real *, real *> idofsAddressRegistry{};
};


class PlasticityRecorder : public AbstractRecorder {
public:
  void setUpContext(LTS &handler, Layer &layer, kernels::LocalData::Loader &loader) {
    currentLoader = &loader;
    AbstractRecorder::setUpContext(handler, layer);
  }

  void record(LTS &handler, Layer &layer) override;
  kernels::LocalData::Loader *currentLoader{nullptr};
};

} // namespace recording
} // namespace initializers
} // namespace seissol

#endif // SEISSOL_RECORDERS_H