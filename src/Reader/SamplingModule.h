// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SAMPLINGMODULE_H_
#define SEISSOL_SRC_READER_SAMPLINGMODULE_H_

#include "Modules/Module.h"

#include <Memory/MemoryAllocator.h>
#include <Modules/Modules.h>
#include <optional>
#include <vector>
namespace seissol::reader {

/**
 * A SeisSol module to sample input at a given interval.
 */
class SamplingModule : public Module {
  public:
  void setName(const std::string& name);

  void setDataSource(std::size_t dataPerSlice,
                     const std::function<void(void*, const std::vector<double>&)>& source);

  void setSpace(std::size_t maxSpace);

  void setup();

  void simulationStart(std::optional<double> time) override;
  void syncPoint(double time) override;

  [[nodiscard]] std::size_t timeSlice(double time) const;
  [[nodiscard]] const void* getData(std::size_t slice) const;
  [[nodiscard]] const void* getData(double time) const;

  private:
  void loadData(double time);

  std::string name_;
  double startTime_{};
  double sampleInterval_{};
  std::size_t dataPerSlice_{0};
  std::function<void(void*, const std::vector<double>&)> sampler_;
  memory::MemkindArray<uint8_t> data_{0, memory::Memkind::Standard};
};

} // namespace seissol::reader
#endif // SEISSOL_SRC_READER_SAMPLINGMODULE_H_
