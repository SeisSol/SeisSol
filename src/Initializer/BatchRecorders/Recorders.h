// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_

#include "Common/ConfigRegistry.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Kernels/Interface.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"

#include <array>
#include <cstddef>
#include <map>
#include <unordered_map>
#include <utility>
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

/// Records the batches of the local integration of a layer of the configuration `Cfg`.
template <typename Cfg>
class LocalIntegrationRecorder : public AbstractRecorder<LTS::LTSVarmap> {
  public:
  using real = Real<Cfg>;

  explicit LocalIntegrationRecorder(double g) : g_(g) {}

  void record(LTS::Layer& layer) override;

  protected:
  double g_{9.81};

  void setUpContext(LTS::Layer& layer) override {
    integratedDofsAddressCounter_ = 0;
    derivativesAddressCounter_ = 0;
    AbstractRecorder::setUpContext(layer);
  }

  private:
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

/// Records the batches of the neighbor integration of a layer of the configuration `Cfg`.
template <typename Cfg>
class NeighIntegrationRecorder : public AbstractRecorder<LTS::LTSVarmap> {
  public:
  using real = Real<Cfg>;

  void record(LTS::Layer& layer) override;

  protected:
  void setUpContext(LTS::Layer& layer) override {
    integratedDofsAddressCounter_ = 0;
    configBoundaryBatches_.assign(builtConfigCount(), {});
    fromCanonical_ = {};
    converted_ = {};
    fromCoupledCanonical_ = {};
    coupledConverted_ = {};
    normalStress_ = {};
    canonicalRegistry_.clear();
    convertedRegistry_.clear();
    configBoundaryScratchCounter_ = 0;
    AbstractRecorder::setUpContext(layer);
  }

  private:
  void recordDofsTimeEvaluation();
  void recordNeighborFluxIntegrals();

  /// Records the face `face` of the cell `cell`, whose neighbor computes in another configuration
  /// (see kernels::ConfigBoundary): the time integral of the neighbor in its configuration if it
  /// provides derivatives, its conversion into the canonical form, and from there into the
  /// configuration of the layer for the side of the neighbor; each once per neighbor (and side).
  void recordConfigBoundaryFace(std::size_t cell, std::size_t face);
  /// Sets the batches recordConfigBoundaryFace collected.
  void recordConfigBoundaryBatches();
  /// The next `bytes` of the scratch of the faces between configurations.
  void* allocateConfigBoundaryScratch(std::size_t bytes);

  std::unordered_map<real*, real*> idofsAddressRegistry_;
  size_t integratedDofsAddressCounter_{0};

  // per configuration of the neighbors
  struct ConfigBoundaryBatches {
    // the derivatives of the neighbors and their time integrals, with the GTS and the LTS relation
    std::array<std::vector<void*>, 2> derivatives;
    std::array<std::vector<void*>, 2> integrals;
    // the time integrals of the neighbors to convert, and their canonical forms
    std::vector<void*> toCanonical;
    std::vector<double*> canonical;
  };
  std::vector<ConfigBoundaryBatches> configBoundaryBatches_;
  // per side of the neighbors: the canonical forms, and their conversions
  std::array<std::vector<double*>, Cell::NumFaces> fromCanonical_;
  std::array<std::vector<real*>, Cell::NumFaces> converted_;
  // the same for the neighbors of the coupled family, with the weights of the normal stress on the
  // face if the cells are of a fluid
  std::array<std::vector<double*>, Cell::NumFaces> fromCoupledCanonical_;
  std::array<std::vector<real*>, Cell::NumFaces> coupledConverted_;
  std::array<std::vector<real*>, Cell::NumFaces> normalStress_;
  // the canonical form of each neighbor, and its conversion per side
  std::unordered_map<const void*, double*> canonicalRegistry_;
  std::map<std::pair<const void*, std::size_t>, real*> convertedRegistry_;
  std::size_t configBoundaryScratchCounter_{0};
};

/// Records the batches of the plasticity of a layer of the configuration `Cfg`.
template <typename Cfg>
class PlasticityRecorder : public AbstractRecorder<LTS::LTSVarmap> {
  public:
  using real = Real<Cfg>;

  protected:
  void setUpContext(LTS::Layer& layer) override { AbstractRecorder::setUpContext(layer); }

  public:
  void record(LTS::Layer& layer) override;
};

/// Records the batches of a layer of dynamic rupture faces of the configuration `Cfg`.
template <typename Cfg>
class DynamicRuptureRecorder : public AbstractRecorder<DynamicRupture::DynrupVarmap> {
  public:
  using real = Real<Cfg>;

  void record(DynamicRupture::Layer& layer) override;

  protected:
  void setUpContext(DynamicRupture::Layer& layer) override {
    AbstractRecorder::setUpContext(layer);
  }

  private:
  void recordSpaceInterpolation();
  std::unordered_map<real*, real*> idofsAddressRegistry_;
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_RECORDERS_H_
