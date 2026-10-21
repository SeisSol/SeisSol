// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "BufferRegistry.h"

#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Writer.h"

#include <algorithm>
#include <cstddef>
#include <unordered_set>
#include <utils/logger.h>
#include <vector>

namespace seissol::io::writer::module {

std::vector<int> BufferRegistry::assign(Writer& writer) {
  if (allocator_.collective()) {
    return assignCollective(writer);
  }

  std::unordered_set<DataSource*> handled;
  std::unordered_set<int> used;
  std::vector<int> ids;

  for (const auto& instruction : writer.getInstructions()) {
    for (auto& dataSource : instruction->dataSources()) {
      if (!dataSource->managed() || handled.find(dataSource.get()) != handled.end()) {
        continue;
      }
      const auto id = idFor(dataSource.get(), used);
      if (id < 0) {
        logError() << "Internal buffer error.";
      }
      used.emplace(id);
      ids.push_back(id);
      dataSource->assignId(id);
      handled.emplace(dataSource.get());
    }
  }

  return ids;
}

std::vector<int> BufferRegistry::assignCollective(Writer& writer) {
  // the sources of the plan, each once, in plan order
  std::vector<DataSource*> sources;
  std::unordered_set<DataSource*> handled;
  for (const auto& instruction : writer.getInstructions()) {
    for (auto& dataSource : instruction->dataSources()) {
      if (dataSource->managed() && handled.emplace(dataSource.get()).second) {
        sources.push_back(dataSource.get());
      }
    }
  }

  // what every source asks for, and the buffer it takes
  struct Use {
    std::vector<BufferSlot>* slots;
    std::size_t slot;
    const void* pointer;
    std::size_t size;
    AdhocBuffer* managed;
  };
  std::vector<Use> uses;
  uses.reserve(sources.size());
  std::size_t passThrough = 0;
  std::size_t managed = 0;
  for (auto* source : sources) {
    if (auto* buffer = dynamic_cast<WriteBuffer*>(source); buffer != nullptr) {
      uses.push_back({&passThroughSlots_,
                      passThrough++,
                      buffer->getLocalPointer(),
                      buffer->getLocalSize(),
                      nullptr});
    } else if (auto* adhoc = dynamic_cast<AdhocBuffer*>(source); adhoc != nullptr) {
      uses.push_back({&managedSlots_, managed++, nullptr, adhoc->getTargetSize(), adhoc});
    } else {
      logError() << "Unsupported buffer type.";
    }
  }
  passThroughSlots_.resize(std::max(passThroughSlots_.size(), passThrough));
  managedSlots_.resize(std::max(managedSlots_.size(), managed));

  // adding and resizing are collective: a buffer that changes on one rank changes on all
  std::vector<int> changed(uses.size());
  for (std::size_t i = 0; i < uses.size(); ++i) {
    const auto& slot = (*uses[i].slots)[uses[i].slot];
    changed[i] =
        (slot.id < 0 || slot.pointer != uses[i].pointer || slot.size != uses[i].size) ? 1 : 0;
  }
  allocator_.agreeOnAny(changed);

  std::vector<int> ids;
  ids.reserve(uses.size());
  for (std::size_t i = 0; i < uses.size(); ++i) {
    const auto& use = uses[i];
    auto& slot = (*use.slots)[use.slot];
    if (slot.id < 0) {
      slot.id = allocator_.addBuffer(use.pointer, use.size);
    } else if (changed[i] != 0) {
      allocator_.resizeBuffer(slot.id, use.pointer, use.size);
    }
    slot.pointer = use.pointer;
    slot.size = use.size;

    if (use.managed != nullptr) {
      // avoid assert in ASYNC by checking the size explicitly
      use.managed->setData(use.size == 0 ? nullptr : allocator_.managedBuffer(slot.id));
    }
    sources[i]->assignId(slot.id);
    ids.push_back(slot.id);
  }

  return ids;
}

int BufferRegistry::idFor(DataSource* source, const std::unordered_set<int>& used) {
  if (auto* passThrough = dynamic_cast<WriteBuffer*>(source); passThrough != nullptr) {
    return idForPassThrough(passThrough, used);
  }
  if (auto* managed = dynamic_cast<AdhocBuffer*>(source); managed != nullptr) {
    return idForManaged(managed, used);
  }
  logError() << "Unsupported buffer type.";
  return -1;
}

int BufferRegistry::idForPassThrough(WriteBuffer* source, const std::unordered_set<int>& used) {
  const auto* pointer = source->getLocalPointer();
  const auto size = source->getLocalSize();

  const auto entry = pointerMap_.find(pointer);
  if (entry == pointerMap_.end()) {
    const BufferPointer added{size, allocator_.addBuffer(pointer, size)};
    pointerMap_[pointer] = added;
    return added.id;
  }

  auto& existing = entry->second;
  if (existing.size != size) {
    if (used.find(existing.id) != used.end()) {
      // it is ok to request the same buffer multiple times, but not with different sizes (that's
      // currently still unsupported)
      logError() << "The same buffer is requested with different sizes. This is not supported at "
                    "the moment.";
    }
    allocator_.resizeBuffer(existing.id, pointer, size);
    existing.size = size;
  }
  return existing.id;
}

int BufferRegistry::idForManaged(AdhocBuffer* source, const std::unordered_set<int>& used) {
  const auto targetSize = source->getTargetSize();

  // recycle an allocation of the same size, unless this write already claimed it. This does of
  // course assume that we don't change the size too often...
  int id = -1;
  for (const auto candidate : bufferMap_[targetSize]) {
    if (used.find(candidate) == used.end()) {
      id = candidate;
      break;
    }
  }
  if (id < 0) {
    id = allocator_.addBuffer(nullptr, targetSize);
    bufferMap_[targetSize].push_back(id);
  }

  // avoid assert in ASYNC by checking targetSize == 0 explicitly
  source->setData(targetSize == 0 ? nullptr : allocator_.managedBuffer(id));
  return id;
}

} // namespace seissol::io::writer::module
