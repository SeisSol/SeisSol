// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "BufferRegistry.h"

#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Writer.h"

#include <unordered_set>
#include <utils/logger.h>
#include <vector>

namespace seissol::io::writer::module {

std::vector<int> BufferRegistry::assign(Writer& writer) {
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
