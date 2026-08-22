// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SamplingModule.h"

#include "Modules/Module.h"

#include <Memory/MemoryAllocator.h>
#include <Modules/Modules.h>
#include <optional>
#include <vector>
namespace seissol::reader {

void SamplingModule::setName(const std::string& name) { this->name_ = name; }

void SamplingModule::setDataSource(
    std::size_t dataPerSlice,
    const std::function<void(void*, const std::vector<double>&)>& source) {
  this->sampler_ = source;
  this->dataPerSlice_ = dataPerSlice;
}

void SamplingModule::setSpace(std::size_t maxSpace) {
  const auto slices = maxSpace / dataPerSlice_;
  data_.resize(slices * dataPerSlice_);
}

void SamplingModule::setup() {
  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);

  const auto slices = data_.size() / dataPerSlice_;
  const auto advance = slices * sampleInterval_;
  setSyncInterval(advance);
}

void SamplingModule::simulationStart(std::optional<double> time) { loadData(time.value_or(0)); }
void SamplingModule::syncPoint(double time) { loadData(time); }

std::size_t SamplingModule::timeSlice(double time) const {
  return static_cast<std::size_t>(std::round((time - startTime_) / sampleInterval_));
}

const void* SamplingModule::getData(std::size_t slice) const {
  return data_.data() + slice * dataPerSlice_;
}

const void* SamplingModule::getData(double time) const { return getData(timeSlice(time)); }

void SamplingModule::loadData(double time) {
  const auto sliceCount = data_.size() / dataPerSlice_;
  std::vector<double> sampleTimes(sliceCount);
  for (std::size_t i = 0; i < sampleTimes.size(); ++i) {
    sampleTimes[i] = time + i * sampleInterval_;
  }
  const auto next = time + sampleTimes.size() * sampleInterval_;

  logInfo() << "Sampling data for" << name_ << "in [" << time << "," << next << ") ...";

  sampler_(data_.data(), sampleTimes);

  logInfo() << "Sampling data for" << name_ << ": Done.";

  startTime_ = time;
}

} // namespace seissol::reader
