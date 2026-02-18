// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_PICKPOINTWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_PICKPOINTWRITER_H_

#include "Modules/Module.h"

#include <functional>

namespace seissol::writer {

class PickpointWriter : public seissol::Module {
  public:
  void syncPoint(double /*currentTime*/) override;
  void simulationStart(std::optional<double> checkpointTime) override;
  void simulationEnd() override;

  void enable(double interval);
  void setupWriter(const std::function<void()>& writer);

  private:
  bool enabled_{false};
  std::function<void()> writer_;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_PICKPOINTWRITER_H_
