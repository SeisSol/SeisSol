#ifndef SEISSOL_RECORDERS_H
#define SEISSOL_RECORDERS_H

#include "DataTypes/ConditionalTable.hpp"
#include "utils/logger.h"
#include <Initializer/LTS.h>
#include <Initializer/DynamicRupture.h>
#include <Initializer/tree/Layer.hpp>
#include "Kernels/common.hpp"
#include <vector>

namespace seissol::initializers::recording {
template <typename LtsT>
class AbstractRecorder {
  public:
  virtual ~AbstractRecorder() = default;

  virtual void record(LtsT& handler, Layer& layer) = 0;

  protected:
  void checkKey(const ConditionalKey& key) {
    if (currentTable->find(key) != currentTable->end()) {
      logError()
          << "Table key conflict detected. Problems with hashing in batch recording subsystem";
    }
  }

  void setUpContext(LtsT& handler, Layer& layer) {
    currentTable = &(layer.getConditionalTable<inner_keys::Wp>());
    currentDrTable = &(layer.getConditionalTable<inner_keys::Dr>());
    currentMaterialTable = &(layer.getConditionalTable<inner_keys::Material>());
    currentIndicesTable = &(layer.getConditionalTable<inner_keys::Indices>());
    currentHandler = &(handler);
    currentLayer = &(layer);
  }

  ConditionalPointersToRealsTable* currentTable{nullptr};
  DrConditionalPointersToRealsTable* currentDrTable{nullptr};
  ConditionalMaterialTable* currentMaterialTable{nullptr};
  ConditionalIndicesTable* currentIndicesTable{nullptr};
  LtsT* currentHandler{nullptr};
  Layer* currentLayer{nullptr};
};

template <typename LtsT>
class CompositeRecorder : public AbstractRecorder<LtsT> {
  public:
  ~CompositeRecorder() override {
    for (auto recorder : concreteRecorders)
      delete recorder;
  }

  void record(LtsT& handler, Layer& layer) override {
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

template <typename Config>
class LocalIntegrationRecorder : public AbstractRecorder<seissol::initializers::LTS<Config>> {
  public:
  using RealT = typename Config::RealT;
  void record(LTS<Config>& handler, Layer& layer) override;

  private:
  void setUpContext(LTS<Config>& handler,
                    Layer& layer,
                    typename kernels::LocalData<Config>::Loader& loader) {
    currentLoader = &loader;
    integratedDofsAddressCounter = 0;
    derivativesAddressCounter = 0;
    AbstractRecorder<LTS<Config>>::setUpContext(handler, layer);
  }

  typename kernels::LocalData<Config>::Loader* currentLoader{nullptr};
  void recordTimeAndVolumeIntegrals();
  void recordFreeSurfaceGravityBc();
  void recordDirichletBc();
  void recordAnalyticalBc();
  void recordLocalFluxIntegral();
  void recordDisplacements();

  std::unordered_map<size_t, RealT*> idofsAddressRegistry{};
  std::vector<RealT*> dQPtrs{};

  size_t integratedDofsAddressCounter{0};
  size_t derivativesAddressCounter{0};
};

template <typename Config>
class NeighIntegrationRecorder : public AbstractRecorder<seissol::initializers::LTS<Config>> {
  public:
  using RealT = typename Config::RealT;
  void record(LTS<Config>& handler, Layer& layer) override;

  private:
  void setUpContext(LTS<Config>& handler,
                    Layer& layer,
                    typename kernels::NeighborData<Config>::Loader& loader) {
    currentLoader = &loader;
    integratedDofsAddressCounter = 0;
    AbstractRecorder<LTS<Config>>::setUpContext(handler, layer);
  }
  void recordDofsTimeEvaluation();
  void recordNeighbourFluxIntegrals();
  typename kernels::NeighborData<Config>::Loader* currentLoader{nullptr};
  std::unordered_map<RealT*, RealT*> idofsAddressRegistry{};
  size_t integratedDofsAddressCounter{0};
};

template <typename Config>
class PlasticityRecorder : public AbstractRecorder<seissol::initializers::LTS<Config>> {
  public:
  using RealT = typename Config::RealT;
  void setUpContext(LTS<Config>& handler,
                    Layer& layer,
                    typename kernels::LocalData<Config>::Loader& loader) {
    currentLoader = &loader;
    AbstractRecorder<LTS<Config>>::setUpContext(handler, layer);
  }

  void record(LTS<Config>& handler, Layer& layer) override;
  typename kernels::LocalData<Config>::Loader* currentLoader{nullptr};
};

template <typename Config>
class DynamicRuptureRecorder
    : public AbstractRecorder<seissol::initializers::DynamicRupture<Config>> {
  public:
  using RealT = typename Config::RealT;
  void record(DynamicRupture<Config>& handler, Layer& layer) override;

  private:
  void setUpContext(DynamicRupture<Config>& handler, Layer& layer) {
    AbstractRecorder<LTS<Config>>::setUpContext(handler, layer);
  }
  void recordDofsTimeEvaluation();
  void recordSpaceInterpolation();
  std::unordered_map<RealT*, RealT*> idofsAddressRegistry{};
};

} // namespace seissol::initializers::recording

#endif // SEISSOL_RECORDERS_H
