// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_KERNELS_BATCHRECORDERS_RECORDERS_H_
#define SEISSOL_SRC_CORE_KERNELS_BATCHRECORDERS_RECORDERS_H_

#include "Core/Kernels/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Core/Kernels/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Core/Memory/Tree/Layer.h"

#include <cstddef>
#include <memory>
#include <utils/logger.h>
#include <vector>

namespace seissol::recording {

class AbstractRecorder {
  public:
  virtual ~AbstractRecorder() = default;

  virtual void record(initializer::Layer& layer) = 0;

  protected:
  void checkKey(const ConditionalKey& key) {
    if (currentTable_->contains(key)) {
      logError()
          << "Table key conflict detected. Problems with hashing in batch recording subsystem";
    }
  }

  virtual void setUpContext(initializer::Layer& layer) {
    currentTable_ = &layer.getConditionalTable();
    currentLayer_ = &layer;
  }

  ConditionalTable* currentTable_{nullptr};
  initializer::Layer* currentLayer_{nullptr};
};

class CompositeRecorder : public AbstractRecorder {
  public:
  void record(initializer::Layer& layer) override {
    for (const auto& recorder : concreteRecorders_) {
      recorder->record(layer);
    }
  }

  void addRecorder(AbstractRecorder* recorder) {
    concreteRecorders_.push_back(std::shared_ptr<AbstractRecorder>(recorder));
  }

  void removeRecorder(size_t recorderIndex) {
    if (recorderIndex < concreteRecorders_.size()) {
      concreteRecorders_.erase(concreteRecorders_.begin() + recorderIndex);
    }
  }

  private:
  std::vector<std::shared_ptr<AbstractRecorder>> concreteRecorders_;
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_CORE_KERNELS_BATCHRECORDERS_RECORDERS_H_
