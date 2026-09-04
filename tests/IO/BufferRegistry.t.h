// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Inference.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Module/BufferRegistry.h"
#include "IO/Writer/Module/WriterModule.h"
#include "IO/Writer/Writer.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {
using namespace seissol::io::writer;
using namespace seissol::io::writer::module;

//! Records what the registry asks of ASYNC, and hands out ascending ids.
class RecordingAllocator : public BufferAllocator {
  public:
  struct Call {
    int id;
    const void* pointer;
    std::size_t size;
    bool resize;
  };

  int addBuffer(const void* pointer, std::size_t size) override {
    const auto id = static_cast<int>(storage_.size());
    storage_.emplace_back(size);
    calls.push_back({id, pointer, size, false});
    return id;
  }

  void resizeBuffer(int id, const void* pointer, std::size_t size) override {
    storage_[id].resize(size);
    calls.push_back({id, pointer, size, true});
  }

  void* managedBuffer(int id) override { return storage_[id].data(); }

  [[nodiscard]] std::size_t bufferCount() const { return storage_.size(); }
  [[nodiscard]] std::size_t bufferSize(int id) const { return storage_[id].size(); }

  std::vector<Call> calls;

  private:
  std::vector<std::vector<char>> storage_;
};

//! A plan with one instruction per data source, in the given order.
Writer planOf(const std::vector<std::shared_ptr<DataSource>>& sources) {
  Writer writer;
  for (std::size_t i = 0; i < sources.size(); ++i) {
    writer.addInstruction(std::make_shared<instructions::Hdf5DataWrite>(
        instructions::Hdf5Location("file.h5", {"group"}),
        "data" + std::to_string(i),
        sources[i],
        seissol::io::datatype::inferDatatype<char>()));
  }
  return writer;
}

std::shared_ptr<DataSource> managedOf(std::size_t count) {
  return GeneratedBuffer::createElementwise<char>(
      count, 1, {}, [](char* target, std::size_t index) { target[0] = static_cast<char>(index); });
}

} // namespace

TEST_CASE("IO/BufferRegistry: sources of equal size do not share an id within a write" *
          doctest::test_suite("io")) {
  RecordingAllocator allocator;
  BufferRegistry registry(allocator);

  // recycling is keyed by size, so two equally sized sources are exactly the case that must not
  // collapse onto one buffer -- they are written in the same step
  auto first = managedOf(64);
  auto second = managedOf(64);
  auto plan = planOf({first, second});

  const auto ids = registry.assign(plan);
  REQUIRE(ids.size() == 2);
  CHECK(ids[0] != ids[1]);
  CHECK(allocator.bufferCount() == 2);

  // and each of them was pointed at its own memory
  CHECK(allocator.managedBuffer(ids[0]) != allocator.managedBuffer(ids[1]));
}

TEST_CASE("IO/BufferRegistry: buffers are recycled across writes" * doctest::test_suite("io")) {
  RecordingAllocator allocator;
  BufferRegistry registry(allocator);

  std::vector<int> seen;
  for (int step = 0; step < 3; ++step) {
    auto plan = planOf({managedOf(64), managedOf(64)});
    const auto ids = registry.assign(plan);
    REQUIRE(ids.size() == 2);
    seen.insert(seen.end(), ids.begin(), ids.end());
  }

  // two allocations, reused in every step
  CHECK(allocator.bufferCount() == 2);
  for (std::size_t i = 2; i < seen.size(); ++i) {
    CHECK(seen[i] == seen[i % 2]);
  }

  // a size that has not been seen before needs a new one
  auto plan = planOf({managedOf(128)});
  const auto ids = registry.assign(plan);
  REQUIRE(ids.size() == 1);
  CHECK(allocator.bufferCount() == 3);
  CHECK(allocator.bufferSize(ids[0]) == 128);
}

TEST_CASE("IO/BufferRegistry: a pass-through buffer keeps the id of its pointer" *
          doctest::test_suite("io")) {
  RecordingAllocator allocator;
  BufferRegistry registry(allocator);

  std::vector<char> data(32);

  auto plan = planOf({WriteBuffer::create(data.data(), data.size())});
  const auto first = registry.assign(plan);
  REQUIRE(first.size() == 1);
  CHECK(allocator.bufferCount() == 1);
  CHECK(allocator.calls.back().pointer == data.data());

  // the same pointer at the same size in a later write must not allocate again
  auto again = planOf({WriteBuffer::create(data.data(), data.size())});
  const auto second = registry.assign(again);
  CHECK(second == first);
  CHECK(allocator.bufferCount() == 1);
  CHECK(allocator.calls.size() == 1);

  // a changed size keeps the id, but has to reach the allocator
  auto grown = planOf({WriteBuffer::create(data.data(), data.size() * 2)});
  const auto third = registry.assign(grown);
  CHECK(third == first);
  CHECK(allocator.bufferCount() == 1);
  REQUIRE(allocator.calls.size() == 2);
  CHECK(allocator.calls.back().resize);
  CHECK(allocator.calls.back().size == data.size() * 2);
}

TEST_CASE("IO/BufferRegistry: one source used twice gets one id" * doctest::test_suite("io")) {
  RecordingAllocator allocator;
  BufferRegistry registry(allocator);

  // the geometry of a mesh output is referenced by several instructions
  auto shared = managedOf(64);
  auto plan = planOf({shared, shared, managedOf(64)});

  const auto ids = registry.assign(plan);
  REQUIRE(ids.size() == 2);
  CHECK(ids[0] != ids[1]);
  CHECK(allocator.bufferCount() == 2);
}

TEST_CASE("IO/BufferRegistry: inline data needs no buffer" * doctest::test_suite("io")) {
  RecordingAllocator allocator;
  BufferRegistry registry(allocator);

  auto plan = planOf({WriteInline::createArray<std::int64_t>({2}, {1, 2}), managedOf(64)});

  const auto ids = registry.assign(plan);
  CHECK(ids.size() == 1);
  CHECK(allocator.bufferCount() == 1);
}

TEST_CASE("IO/BufferRegistry: the ids reach the plan" * doctest::test_suite("io")) {
  RecordingAllocator allocator;
  BufferRegistry registry(allocator);

  auto plan = planOf({managedOf(64), managedOf(32)});
  const auto ids = registry.assign(plan);
  REQUIRE(ids.size() == 2);

  // the executor finds the data by the id in the serialised plan, so the two have to agree
  const auto serialized = plan.serialize();
  for (const auto id : ids) {
    CHECK(serialized.find("id: " + std::to_string(id)) != std::string::npos);
  }
}

} // namespace seissol::unit_test

// ---------------------------------------------------------------------------
// Resuming the output numbering after a restart
// ---------------------------------------------------------------------------

TEST_CASE("IO/WriterModule: the output numbering resumes after a restart" *
          doctest::test_suite("io")) {
  using seissol::io::writer::module::outputCountBefore;

  // a fresh run has nothing behind it
  CHECK(outputCountBefore(0.0, 0.5) == 0);

  // outputs happen at 0, 0.5, 1.0, ...; a checkpoint at 0.5 has two of them behind it
  CHECK(outputCountBefore(0.5, 0.5) == 2);
  CHECK(outputCountBefore(1.0, 0.5) == 3);
  CHECK(outputCountBefore(2.0, 0.5) == 5);

  // a checkpoint between two outputs belongs to the earlier one
  CHECK(outputCountBefore(0.7, 0.5) == 2);
  CHECK(outputCountBefore(0.99, 0.5) == 2);

  // a time that should sit exactly on an output can arrive a rounding error below it, and must
  // not be counted as one output short -- that would overwrite the last file of the previous run
  double accumulated = 0;
  for (int step = 0; step < 30; ++step) {
    accumulated += 0.1;
  }
  CHECK(accumulated != 3.0);
  CHECK(outputCountBefore(accumulated, 0.1) == 31);

  // nothing sensible to resume from
  CHECK(outputCountBefore(3.0, 0.0) == 0);
  CHECK(outputCountBefore(-1.0, 0.5) == 0);
}
