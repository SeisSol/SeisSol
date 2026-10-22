// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Solver/TimeStepping/Plan/SuperStepRecorder.h"

#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::unit_test {

#ifdef ACL_DEVICE

namespace supersteprecorder {

/// The device work of a super-timestep: copies `source` to `target` on the stream.
class CopyOnDevice {
  public:
  CopyOnDevice()
      : api_(::device::DeviceInstance::instance().api()), stream_(nullptr), source_(nullptr),
        target_(nullptr) {
    stream_ = api_.createStream();
    source_ = static_cast<int*>(api_.allocGlobMem(sizeof(int)));
    target_ = static_cast<int*>(api_.allocGlobMem(sizeof(int)));
    setSource(0);
    const int zero = 0;
    api_.copyTo(target_, &zero, sizeof(int));
  }
  ~CopyOnDevice() {
    api_.syncDevice();
    api_.freeGlobMem(source_);
    api_.freeGlobMem(target_);
    api_.destroyGenericStream(stream_);
  }
  CopyOnDevice(const CopyOnDevice&) = delete;
  CopyOnDevice(CopyOnDevice&&) = delete;
  CopyOnDevice& operator=(const CopyOnDevice&) = delete;
  CopyOnDevice& operator=(CopyOnDevice&&) = delete;

  void setSource(int value) { api_.copyTo(source_, &value, sizeof(int)); }

  void enqueue() { api_.copyBetweenAsync(target_, source_, sizeof(int), stream_); }

  /// the target, once the device is past everything
  [[nodiscard]] int target() {
    api_.syncDevice();
    int value = -1;
    api_.copyFrom(&value, target_, sizeof(int));
    return value;
  }

  [[nodiscard]] std::vector<void*> streams() const { return {stream_}; }

  private:
  ::device::AbstractAPI& api_;
  void* stream_{nullptr};
  int* source_{nullptr};
  int* target_{nullptr};
};

} // namespace supersteprecorder

TEST_CASE("A recorded super-timestep replays its work, on the data of the time of the replay" *
          doctest::test_suite("solver")) {
  using namespace supersteprecorder;
  auto& device = ::device::DeviceInstance::instance();
  device.api().setDevice(0);
  device.api().initialize();
  if (!solver::SuperStepRecorder::available()) {
    MESSAGE("this device cannot record graphs; not tested");
    return;
  }

  CopyOnDevice work;
  solver::SuperStepRecorder recorder;
  const solver::SuperStepRecorder::Key key{0.5, 1.0};
  const solver::SuperStepRecorder::Key otherKey{0.25, 1.0};
  CHECK_FALSE(recorder.has(key));

  // a recording only runs once replayed
  work.setSource(1);
  recorder.beginRecording({}, work.streams());
  work.enqueue();
  recorder.endRecording(key, {}, work.streams());
  CHECK(recorder.has(key));
  CHECK_FALSE(recorder.has(otherKey));
  CHECK(work.target() == 0);
  recorder.replay(key, {}, work.streams());
  CHECK(work.target() == 1);

  // a replay takes the data as it is then
  work.setSource(2);
  recorder.beginReplay({}, work.streams());
  recorder.replay(key, {}, work.streams());
  CHECK(work.target() == 2);

  // launched once, a recording runs, but is not kept
  work.setSource(3);
  recorder.beginRecording({}, work.streams());
  work.enqueue();
  recorder.endRecordingAndLaunch({}, work.streams());
  CHECK(work.target() == 3);
  CHECK_FALSE(recorder.has(otherKey));

  // an epoch that has replayed a recording keeps it; the next one without a replay lets it go
  recorder.forgetUnused();
  CHECK(recorder.has(key));
  recorder.forgetUnused();
  CHECK_FALSE(recorder.has(key));

  recorder.dispose();
}

#endif

} // namespace seissol::unit_test
