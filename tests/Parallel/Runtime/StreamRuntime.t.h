// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/Runtime/Event.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cstddef>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::unit_test {

#ifdef ACL_DEVICE

namespace streamruntime {

/// a few numbers on the device
class DeviceNumbers {
  public:
  static constexpr std::size_t Count = 8;

  DeviceNumbers()
      : api_(::device::DeviceInstance::instance().api()),
        data_(static_cast<int*>(api_.allocGlobMem(Count * sizeof(int)))) {
    set({});
  }
  ~DeviceNumbers() {
    api_.syncDevice();
    api_.freeGlobMem(data_);
  }
  DeviceNumbers(const DeviceNumbers&) = delete;
  DeviceNumbers(DeviceNumbers&&) = delete;
  DeviceNumbers& operator=(const DeviceNumbers&) = delete;
  DeviceNumbers& operator=(DeviceNumbers&&) = delete;

  void set(const std::array<int, Count>& values) {
    api_.syncDevice();
    api_.copyTo(data_, values.data(), Count * sizeof(int));
  }

  [[nodiscard]] std::array<int, Count> get() {
    api_.syncDevice();
    std::array<int, Count> values{};
    api_.copyFrom(values.data(), data_, Count * sizeof(int));
    return values;
  }

  /// copies one number to another one on the stream
  void copy(std::size_t to, std::size_t from, void* stream) {
    api_.copyBetweenAsync(data_ + to, data_ + from, sizeof(int), stream);
  }

  private:
  ::device::AbstractAPI& api_;
  int* data_{nullptr};
};

} // namespace streamruntime

TEST_CASE("Stream runtimes that build one graph wait for each other through its edges" *
          doctest::test_suite("parallel")) {
  using namespace streamruntime;
  const auto& device = ::device::DeviceInstance::instance();
  device.api().setDevice(0);
  device.api().initialize();
  if (!device.api().isCapableOfGraphNodes()) {
    MESSAGE("this device cannot build graphs node by node; not tested");
    return;
  }

  // the branches of envMany as nodes of their own, and one after another
  for (const bool branches : {true, false}) {
    CAPTURE(branches);
    DeviceNumbers numbers;
    parallel::runtime::StreamRuntime root;
    parallel::runtime::StreamRuntime first;
    parallel::runtime::StreamRuntime second(2);

    root.beginGraph();
    first.attach(root);
    second.attach(root, branches);

    // number 0 to 1 on the first runtime, then on to 2 on the second one, which waits for it; and 4
    // to 5 on the first one in between, while the second one has a node open
    numbers.copy(1, 0, first.stream());
    second.eventSync(first.eventRecord());
    numbers.copy(2, 1, second.stream());
    numbers.copy(5, 4, first.stream());
    // independent branches of the second runtime take 2 on to 3 and 6, and come back together
    second.envMany(
        2, [&](void* stream, std::size_t branch) { numbers.copy(branch == 0 ? 3 : 6, 2, stream); });
    numbers.copy(7, 3, second.stream());
    const auto graph = root.endGraph();

    // nothing has run yet, and the runtimes run their work directly again
    using Numbers = std::array<int, DeviceNumbers::Count>;
    numbers.set({1, 0, 0, 0, 2, 0, 0, 0});
    CHECK(numbers.get() == Numbers{1, 0, 0, 0, 2, 0, 0, 0});

    root.launchGraph(graph);
    CHECK(numbers.get() == Numbers{1, 1, 1, 1, 2, 2, 1, 1});

    // a replay works on the data of the time of the replay
    numbers.set({3, 0, 0, 0, 4, 0, 0, 0});
    root.launchGraph(graph);
    CHECK(numbers.get() == Numbers{3, 3, 3, 3, 4, 4, 3, 3});

    numbers.copy(0, 5, first.stream());
    CHECK(numbers.get()[0] == 4);

    root.dispose();
    first.dispose();
    second.dispose();
  }
}

TEST_CASE("A stream runtime enqueues a deferred wait with its next work" *
          doctest::test_suite("parallel")) {
  using namespace streamruntime;
  const auto& device = ::device::DeviceInstance::instance();
  device.api().setDevice(0);
  device.api().initialize();

  DeviceNumbers numbers;
  parallel::runtime::StreamRuntime before;
  parallel::runtime::StreamRuntime after;

  numbers.set({1, 0, 0, 0, 0, 0, 0, 0});
  numbers.copy(1, 0, before.stream());
  {
    const auto event = before.eventRecord();
    after.deferWait(event);
    CHECK(after.defersWaitFor(event));
    CHECK_FALSE(after.defersWaitFor(parallel::runtime::Event()));

    // nothing enqueued, nothing waited for; then the next work waits
    numbers.copy(2, 1, after.stream());
    CHECK_FALSE(after.defersWaitFor(event));
    CHECK(numbers.get()[2] == 1);
  }

  // or right away, on request; either way, the runtime lets go of the event
  {
    const auto event = before.eventRecord();
    after.deferWait(event);
    after.enqueueDeferredWait();
    CHECK_FALSE(after.defersWaitFor(event));
  }
  CHECK(before.eventPool().referenced() == 0);

  before.dispose();
  after.dispose();
}

#endif

} // namespace seissol::unit_test
