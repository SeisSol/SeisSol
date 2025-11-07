// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_CHECKPOINT_CHECKPOINTMANAGER_H_
#define SEISSOL_SRC_IO_INSTANCE_CHECKPOINT_CHECKPOINTMANAGER_H_

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Writer/Instructions/Data.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"

#include <map>
#include <string>
#include <utils/logger.h>

namespace seissol::io::writer {
class Writer;
} // namespace seissol::io::writer

namespace seissol::io::instance::checkpoint {

struct CheckpointVariable {
  std::string name;
  void* data;
  std::shared_ptr<datatype::Datatype> datatype;
  std::shared_ptr<datatype::Datatype> memoryDatatype;
  std::optional<std::function<void(void*, const void*)>> pack;
  std::optional<std::function<void(void*, const void*)>> unpack;
};

struct CheckpointTree {
  std::string name;
  std::size_t cells;
  std::vector<std::size_t> ids;
  std::vector<CheckpointVariable> variables;
};

class CheckpointManager {
  public:
  template <typename VarmapT>
  void registerTree(const std::string& name,
                    initializer::Storage<VarmapT>& storage,
                    const std::vector<std::size_t>& ids) {
    dataRegistry[&storage].name = name;
    dataRegistry[&storage].cells = storage.size(Ghost);
    dataRegistry[&storage].ids = ids;
  }

  template <typename HandleT, typename VarmapT>
  void registerData(const std::string& name,
                    initializer::Storage<VarmapT>& storage,
                    const HandleT& var) {
    if (storage.info(var).mask != initializer::LayerMask(Ghost)) {
      logError() << "Invalid layer mask for a checkpointing variable (i.e.: NYI).";
    }
    dataRegistry[&storage].variables.emplace_back(
        CheckpointVariable{name,
                           storage.var(var),
                           datatype::inferDatatype<typename HandleT::Type>(),
                           datatype::inferDatatype<typename HandleT::Type>(),
                           {},
                           {}});
  }

  template <typename StorageT, typename VarmapT>
  void registerData(const std::string& name, initializer::Storage<VarmapT>& storage) {
    if (storage.template info<StorageT>().mask != initializer::LayerMask(Ghost)) {
      logError() << "Invalid layer mask for a checkpointing variable (i.e.: NYI).";
    }
    dataRegistry[&storage].variables.emplace_back(
        CheckpointVariable{name,
                           storage.template var<StorageT>(),
                           datatype::inferDatatype<typename StorageT::Type>(),
                           datatype::inferDatatype<typename StorageT::Type>(),
                           {},
                           {}});
  }

  template <typename S, typename T, typename VarmapT>
  void registerTransformedData(const std::string& name,
                               initializer::Storage<VarmapT>& storage,
                               initializer::Variable<T> var,
                               const std::function<void(void*, const void*)>& pack,
                               const std::function<void(void*, const void*)>& unpack) {
    if (var.mask != initializer::LayerMask(Ghost)) {
      logError() << "Invalid layer mask for a checkpointing variable (i.e.: NYI).";
    }
    dataRegistry[&storage].variables.emplace_back(CheckpointVariable{name,
                                                                     storage.var(var),
                                                                     datatype::inferDatatype<S>(),
                                                                     datatype::inferDatatype<T>(),
                                                                     pack,
                                                                     unpack});
  }

  template <std::size_t Pad, std::size_t Nopad, typename T, std::size_t Npad, typename VarmapT>
  void registerPaddedData(const std::string& name,
                          initializer::Storage<VarmapT>& storage,
                          initializer::Variable<T[Npad]> var) {
    constexpr std::size_t Lines = (Npad / Pad);
    constexpr std::size_t Nnopad = Lines * Nopad;
    constexpr std::size_t Nmin = Nopad < Pad ? Nopad : Pad;
    using Tpad = T[Npad];
    using Tnopad = T[Nnopad];
    registerTransformedData<Tnopad, Tpad>(
        name,
        storage,
        var,
        [](void* nopadV, const void* padV) {
          auto* nopad = reinterpret_cast<T*>(nopadV);
          const auto* pad = reinterpret_cast<const T*>(padV);
          for (std::size_t i = 0; i < Lines; ++i) {
            std::memcpy(nopad + i * Nopad, pad + i * Pad, Nmin * sizeof(T));
            if constexpr (Pad < Nopad) {
              std::memset(nopad + i * Nopad + Pad, 0, Nopad - Pad);
            }
          }
        },
        [](void* padV, const void* nopadV) {
          const auto* nopad = reinterpret_cast<const T*>(nopadV);
          auto* pad = reinterpret_cast<T*>(padV);
          for (std::size_t i = 0; i < Lines; ++i) {
            std::memcpy(pad + i * Pad, nopad + i * Nopad, Nmin * sizeof(T));
            if constexpr (Pad > Nopad) {
              std::memset(pad + i * Pad + Nopad, 0, Pad - Nopad);
            }
          }
        });
  }

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter();

  double loadCheckpoint(const std::string& file);

  private:
  std::map<void*, CheckpointTree> dataRegistry;
};

} // namespace seissol::io::instance::checkpoint

#endif // SEISSOL_SRC_IO_INSTANCE_CHECKPOINT_CHECKPOINTMANAGER_H_
