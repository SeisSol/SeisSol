// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_

#include "Common/Constants.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Kernels/Interface.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"

#include <array>
#include <cstddef>
#include <memory>
#include <utils/logger.h>
#include <vector>

namespace seissol::recording {

/// The field the faces of `cell` of `layer` lift into on the device: its degrees of freedom, or
/// the extended field of a solver that keeps the memory variables apart, which the layer keeps as a
/// scratchpad.
real* liftedFieldOf(LTS::Layer& layer, std::size_t cell);

/// What they lift into instead where the cell is weight-adjusted: its flux buffer, to which it
/// applies the inverse weight of its medium once they are all in (see WeightAdjustable).
real* fluxTargetOf(LTS::Layer& layer, std::size_t cell);

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

  virtual void setUpContext(initializer::Layer<VarmapT>& layer) {
    currentTable_ = &layer.template getConditionalTable<inner_keys::Wp>();
    currentDrTable_ = &layer.template getConditionalTable<inner_keys::Dr>();
    currentIndicesTable_ = &layer.template getConditionalTable<inner_keys::Indices>();
    currentLayer_ = &layer;
  }

  ConditionalPointersToRealsTable* currentTable_{nullptr};
  DrConditionalPointersToRealsTable* currentDrTable_{nullptr};
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
  explicit LocalIntegrationRecorder(double g) : g_(g) {}

  void record(LTS::Layer& layer) override;

  protected:
  double g_{9.81};

  private:
  void recordTimeAndVolumeIntegrals();
  void recordFreeSurfaceGravityBc();
  void recordDirichletBc();
  void recordAnalyticalBc(LTS::Layer& layer);
  void recordFaceTraces();
  void recordDisplacements();

  // the time integral and the time derivatives of every cell
  std::vector<real*> idofsPtrs_;
  std::vector<real*> dQPtrs_;
  // what the faces of every cell lift into, see fluxTargetOf
  std::vector<real*> fluxTargets_;
};

class NeighIntegrationRecorder : public AbstractRecorder<LTS::LTSVarmap> {
  public:
  void record(LTS::Layer& layer) override;

  private:
  void recordNeighborTraces();
  void recordFaceCorrectors();
  void recordDynamicRuptureFluxes();

  // what the corrector of each face of every cell reads as the neighbor's trace
  std::vector<std::array<real*, Cell::NumFaces>> neighborTraces_;
};

class PlasticityRecorder : public AbstractRecorder<LTS::LTSVarmap> {
  protected:
  void setUpContext(LTS::Layer& layer) override { AbstractRecorder::setUpContext(layer); }

  public:
  void record(LTS::Layer& layer) override;
};

class DynamicRuptureRecorder : public AbstractRecorder<DynamicRupture::DynrupVarmap> {
  public:
  void record(DynamicRupture::Layer& layer) override;

  private:
  void recordSpaceInterpolation();
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
