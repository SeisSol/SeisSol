#ifndef SEISSOL_RECORDERS_H
#define SEISSOL_RECORDERS_H

#include "DataTypes/ConditionalTable.hpp"
#include "utils/logger.h"
#include <Initializer/LTS.h>
#include <Initializer/DynamicRupture.h>
#include <Initializer/tree/Layer.hpp>
#include <Kernels/Interface.hpp>
#include <vector>

namespace seissol {
namespace initializers {
namespace recording {


template<typename LtsT>
class AbstractRecorder {
public:
  virtual ~AbstractRecorder() = default;

  virtual void record(LtsT &handler, Layer &layer) = 0;

protected:
  void checkKey(const ConditionalKey &key) {
    if (currentTable->find(key) != currentTable->end()) {
      logError()
          << "Table key conflict detected. Problems with hashing in batch recording subsystem";
    }
  }

  void setUpContext(LtsT &handler, Layer &layer) {
    currentTable = &(layer.getCondBatchTable());
    currentHandler = &(handler);
    currentLayer = &(layer);
  }

  ConditionalBatchTableT *currentTable{nullptr};
  LtsT *currentHandler{nullptr};
  Layer *currentLayer{nullptr};
};

template<typename LtsT>
class CompositeRecorder : public AbstractRecorder<LtsT> {
public:
  ~CompositeRecorder() override {
    for (auto recorder : concreteRecorders)
      delete recorder;
  }

  void record(LtsT &handler, Layer &layer) override {
    for (auto recorder : concreteRecorders) {
      recorder->record(handler, layer);
    }
  }

  void addRecorder(AbstractRecorder<LtsT> *recorder) {
    concreteRecorders.push_back(recorder);
  }

  void removeRecorder(size_t recorderIndex) {
    if (recorderIndex < concreteRecorders.size()) {
      concreteRecorders.erase(concreteRecorders.begin() + recorderIndex);
    }
  }

private:
  std::vector<AbstractRecorder<LtsT> *> concreteRecorders{};
};


class LocalIntegrationRecorder : public AbstractRecorder<seissol::initializers::LTS> {
public:
  void record(LTS &handler, Layer &layer) override;

private:
  void setUpContext(LTS &handler, Layer &layer, kernels::LocalData::Loader &loader) {
    currentLoader = &loader;
    integratedDofsAddressCounter = 0;
    derivativesAddressCounter = 0;
    AbstractRecorder::setUpContext(handler, layer);
  }

  kernels::LocalData::Loader *currentLoader{nullptr};
  void recordTimeAndVolumeIntegrals();
  void recordLocalFluxIntegral();
  void recordDisplacements();
  std::unordered_map<size_t, real *> idofsAddressRegistry{};
  size_t integratedDofsAddressCounter{0};
  size_t derivativesAddressCounter{0};
};


class NeighIntegrationRecorder : public AbstractRecorder<seissol::initializers::LTS> {
public:
  void record(LTS &handler, Layer &layer) override;
private:
  void setUpContext(LTS &handler, Layer &layer, kernels::NeighborData::Loader &loader) {
    currentLoader = &loader;
    integratedDofsAddressCounter = 0;
    AbstractRecorder::setUpContext(handler, layer);
  }
  void recordDofsTimeEvaluation();
  void recordNeighbourFluxIntegrals();
  kernels::NeighborData::Loader *currentLoader{nullptr};
  std::unordered_map<real *, real *> idofsAddressRegistry{};
  size_t integratedDofsAddressCounter{0};
};


class PlasticityRecorder : public AbstractRecorder<seissol::initializers::LTS> {
public:
  void setUpContext(LTS &handler, Layer &layer, kernels::LocalData::Loader &loader) {
    currentLoader = &loader;
    AbstractRecorder::setUpContext(handler, layer);
  }

  void record(LTS &handler, Layer &layer) override;
  kernels::LocalData::Loader *currentLoader{nullptr};
};

class DynamicRuptureRecorder : public AbstractRecorder<seissol::initializers::DynamicRupture> {
public:
  void record(DynamicRupture &handler, Layer &layer) override;
private:
  void setUpContext(DynamicRupture &handler, Layer &layer) {
    AbstractRecorder::setUpContext(handler, layer);
  }
  void recordDofsTimeEvaluation();
  void recordSpaceInterpolation();
  std::unordered_map<real *, real *> idofsAddressRegistry{};
};

} // namespace recording
} // namespace initializers
} // namespace seissol

#endif // SEISSOL_RECORDERS_H
