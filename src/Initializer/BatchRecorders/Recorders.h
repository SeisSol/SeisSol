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
  void checkKey(const ConditionalBatchTableT &table, const ConditionalKey &key) {
    if (table.find(key) != table.end()) {
      logError()
          << "Table key conflict detected. Problems with hashing in batch recording subsystem";
    }
  }

  void loadContext(LTS &handler, Layer &layer, kernels::LocalData::Loader& loader) {
    currentLoader = &loader;
    currentTable = &(layer.getCondBatchTable());
    currentHandler = &(handler);
    currentLayer = &(layer);
  }

  kernels::LocalData::Loader *currentLoader{nullptr};
  ConditionalBatchTableT *currentTable{nullptr};
  LTS *currentHandler{nullptr};
  Layer *currentLayer{nullptr};
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
  void recordTimeIntegral();
  void recordVolumeIntegral();
  void recordLocalFluxIntegral();
  void recordDisplacements();
};


class NeighIntegrationRecorder : public AbstractRecorder {
public:
  void record(LTS &handler, Layer &layer) override;
};


class PlasticityRecorder : public AbstractRecorder {
public:
  void record(LTS &handler, Layer &layer) override;
};

} // namespace recording
} // namespace initializers
} // namespace seissol

#endif // SEISSOL_RECORDERS_H
