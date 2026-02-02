// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_

#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Kernels/Interface.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"

#include <utils/logger.h>
#include <vector>

namespace seissol::recording {

template <typename VarmapT>
class AbstractRecorder {
  public:
  virtual ~AbstractRecorder() = default;

  virtual void record(initializer::Layer<VarmapT>& layer) = 0;

  protected:
  void checkKey(const ConditionalKey& key) {
    if (currentTable_->find(key) != currentTable_->end()) {
      logError()
          << "Table key conflict detected. Problems with hashing in batch recording subsystem";
    }
  }

  void setUpContext(initializer::Layer<VarmapT>& layer) {
    currentTable_ = &(layer.template getConditionalTable<inner_keys::Wp>());
    currentDrTable_ = &(layer.template getConditionalTable<inner_keys::Dr>());
    currentMaterialTable_ = &(layer.template getConditionalTable<inner_keys::Material>());
    currentIndicesTable_ = &(layer.template getConditionalTable<inner_keys::Indices>());
    currentLayer_ = &(layer);
  }

  ConditionalPointersToRealsTable* currentTable_{nullptr};
  DrConditionalPointersToRealsTable* currentDrTable_{nullptr};
  ConditionalMaterialTable* currentMaterialTable_{nullptr};
  ConditionalIndicesTable* currentIndicesTable_{nullptr};
  initializer::Layer<VarmapT>* currentLayer_{nullptr};
};

template <typename VarmapT>
class CompositeRecorder : public AbstractRecorder<VarmapT> {
  public:
  void record(initializer::Layer<VarmapT>& layer) override {
    for (auto& recorder : concreteRecorders_) {
      recorder->record(layer);
    }
  }

  void addRecorder(AbstractRecorder<VarmapT>* recorder) {
    concreteRecorders_.push_back(std::shared_ptr<AbstractRecorder<VarmapT>>(recorder));
  }

  void removeRecorder(size_t recorderIndex) {
    if (recorderIndex < concreteRecorders_.size()) {
      concreteRecorders_.erase(concreteRecorders_.begin() + recorderIndex);
    }
  }

  private:
  std::vector<std::shared_ptr<AbstractRecorder<VarmapT>>> concreteRecorders_;
};

class LocalIntegrationRecorder : public AbstractRecorder<LTS::LTSVarmap> {
  public:
  void record(LTS::Layer& layer) override;

  private:
  void setUpContext(LTS::Layer& layer) {
    integratedDofsAddressCounter_ = 0;
    derivativesAddressCounter_ = 0;
    AbstractRecorder::setUpContext(layer);
  }

  void recordTimeAndVolumeIntegrals();
  void recordFreeSurfaceGravityBc();
  void recordDirichletBc();
  void recordAnalyticalBc(LTS::Layer& layer);
  void recordLocalFluxIntegral();
  void recordDisplacements();

  std::unordered_map<size_t, real*> idofsAddressRegistry_;
  std::vector<real*> dQPtrs_;

  size_t integratedDofsAddressCounter_{0};
  size_t derivativesAddressCounter_{0};
};

class NeighIntegrationRecorder : public AbstractRecorder<LTS::LTSVarmap> {
  public:
  void record(LTS::Layer& layer) override;

  private:
  void setUpContext(LTS::Layer& layer) {
    integratedDofsAddressCounter_ = 0;
    AbstractRecorder::setUpContext(layer);
  }
  void recordDofsTimeEvaluation();
  void recordNeighborFluxIntegrals();
  std::unordered_map<real*, real*> idofsAddressRegistry_;
  size_t integratedDofsAddressCounter_{0};
};

class PlasticityRecorder : public AbstractRecorder<LTS::LTSVarmap> {
  public:
  void setUpContext(LTS::Layer& layer) { AbstractRecorder::setUpContext(layer); }

  void record(LTS::Layer& layer) override;
};

class DynamicRuptureRecorder : public AbstractRecorder<DynamicRupture::DynrupVarmap> {
  public:
  void record(DynamicRupture::Layer& layer) override;

  private:
  void setUpContext(DynamicRupture::Layer& layer) { AbstractRecorder::setUpContext(layer); }
  void recordDofsTimeEvaluation();
  void recordSpaceInterpolation();
  std::unordered_map<real*, real*> idofsAddressRegistry_;
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
