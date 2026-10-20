// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_PICKPOINTWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_PICKPOINTWRITER_H_

#include "Memory/HostData.h"
#include "Modules/Module.h"

#include <functional>

namespace seissol::writer {

class PickpointWriter : public seissol::Module {
  public:
  void captureSyncPoint(double /*currentTime*/) override;
  void processSyncPoint(double /*currentTime*/) override;
  /// The fault receivers take their samples into buffers of their own.
  [[nodiscard]] initializer::HostData hostData() const override { return {}; }
  /// The files get written from what the capture took.
  [[nodiscard]] bool processesAsynchronously() const override { return true; }
  void simulationStart(std::optional<double> checkpointTime) override;
  void simulationEnd() override;

  void enable(double interval);
  /// `capture` takes the samples cached so far, `write` writes what it took.
  void setupWriter(const std::function<void()>& capture, const std::function<void()>& write);

  private:
  bool enabled_{false};
  std::function<void()> capture_;
  std::function<void()> write_;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_PICKPOINTWRITER_H_
