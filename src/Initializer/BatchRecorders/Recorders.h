// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
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

class LocalIntegrationRecorder : public AbstractRecorder<seissol::initializer::LTS> {
  public:
  void record(LTS& handler, Layer& layer) override;

  private:
  void setUpContext(LTS& handler,
                    Layer& layer,
                    kernels::LocalData::Loader& loader,
                    kernels::LocalData::Loader& loaderHost) {
    currentLoader = &loader;
    currentLoaderHost = &loaderHost;
    integratedDofsAddressCounter = 0;
    derivativesAddressCounter = 0;
    AbstractRecorder::setUpContext(handler, layer);
  }

  kernels::LocalData::Loader* currentLoader{nullptr};
  kernels::LocalData::Loader* currentLoaderHost{nullptr};
  void recordTimeAndVolumeIntegrals();
  void recordFreeSurfaceGravityBc();
  void recordDirichletBc();
  void recordAnalyticalBc(LTS& handler, Layer& layer);
  void recordLocalFluxIntegral();
  void recordDisplacements();

  std::unordered_map<size_t, real*> idofsAddressRegistry{};
  std::vector<real*> dQPtrs{};

  size_t integratedDofsAddressCounter{0};
  size_t derivativesAddressCounter{0};
};

class NeighIntegrationRecorder : public AbstractRecorder<seissol::initializer::LTS> {
  public:
  void record(LTS& handler, Layer& layer) override;

  private:
  void setUpContext(LTS& handler,
                    Layer& layer,
                    kernels::NeighborData::Loader& loader,
                    kernels::NeighborData::Loader& loaderHost) {
    currentLoader = &loader;
    currentLoaderHost = &loaderHost;
    integratedDofsAddressCounter = 0;
    AbstractRecorder::setUpContext(handler, layer);
  }
  void recordDofsTimeEvaluation();
  void recordNeighbourFluxIntegrals();
  kernels::NeighborData::Loader* currentLoader{nullptr};
  kernels::NeighborData::Loader* currentLoaderHost{nullptr};
  std::unordered_map<real*, real*> idofsAddressRegistry{};
  size_t integratedDofsAddressCounter{0};
};

class PlasticityRecorder : public AbstractRecorder<seissol::initializer::LTS> {
  public:
  void setUpContext(LTS& handler,
                    Layer& layer,
                    kernels::LocalData::Loader& loader,
                    kernels::LocalData::Loader& loaderHost) {
    currentLoader = &loader;
    currentLoaderHost = &loaderHost;
    AbstractRecorder::setUpContext(handler, layer);
  }

  void record(LTS& handler, Layer& layer) override;
  kernels::LocalData::Loader* currentLoader{nullptr};
  kernels::LocalData::Loader* currentLoaderHost{nullptr};
};

class DynamicRuptureRecorder : public AbstractRecorder<seissol::initializer::DynamicRupture> {
  public:
  void record(DynamicRupture& handler, Layer& layer) override;

  private:
  void setUpContext(DynamicRupture& handler, Layer& layer) {
    AbstractRecorder::setUpContext(handler, layer);
  }
  void recordDofsTimeEvaluation();
  void recordSpaceInterpolation();
  std::unordered_map<real*, real*> idofsAddressRegistry{};
};

} // namespace seissol::initializer::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
