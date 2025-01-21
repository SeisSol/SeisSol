// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_PARALLEL_PIN_H_
#define SEISSOL_SRC_PARALLEL_PIN_H_

#include "Common/IntegerMaskParser.h"
#include "async/as/Pin.h"
#include <deque>
#include <sched.h>
#include <string>

namespace seissol::parallel {

class Pinning {
  private:
  async::as::CpuMask openmpMask{};
  async::as::CpuMask onlineMask{};
  IntegerMaskParser::MaskType parsedFreeCPUsMask;

  public:
  Pinning();

  static std::deque<bool> parseOnlineCpuMask(std::string s, unsigned numberOfConfiguredCpus);
  async::as::CpuMask computeOnlineCpuMask();
  [[nodiscard]] static async::as::CpuMask getWorkerUnionMask();
  [[nodiscard]] async::as::CpuMask getFreeCPUsMask() const;
  static bool freeCPUsMaskEmpty(const async::as::CpuMask& mask);
  [[nodiscard]] async::as::CpuMask getOnlineMask() const;
  [[nodiscard]] static bool areAllCpusOnline();
  void pinToFreeCPUs() const;
  static std::string maskToString(const async::as::CpuMask& mask);
  [[nodiscard]] static async::as::CpuMask getNodeMask();
  void checkEnvVariables();
};

} // namespace seissol::parallel

#endif // SEISSOL_SRC_PARALLEL_PIN_H_
