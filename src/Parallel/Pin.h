// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_PARALLEL_PIN_H_
#define SEISSOL_SRC_PARALLEL_PIN_H_

#include "Common/IntegerMaskParser.h"

#include <async/as/Pin.h>
#include <deque>
#include <optional>
#include <string>
#include <utils/env.h>

namespace seissol::parallel {

class Pinning {
  private:
  async::as::CpuMask openmpMask_{};
  async::as::CpuMask onlineMask_{};
  IntegerMaskParser::MaskType parsedFreeCPUsMask_;

  public:
  Pinning();

  static std::deque<bool> parseOnlineCpuMask(std::string s, unsigned numberOfConfiguredCpus);
  async::as::CpuMask computeOnlineCpuMask();
  [[nodiscard]] static async::as::CpuMask getWorkerUnionMask();
  [[nodiscard]] async::as::CpuMask getFreeCPUsMask() const;
  static bool freeCPUsMaskEmpty(const async::as::CpuMask& mask);
  /**
   * @return The lowest CPU in @p mask, or none if the mask is empty.
   */
  [[nodiscard]] static std::optional<int> firstCpu(const async::as::CpuMask& mask);
  [[nodiscard]] async::as::CpuMask getOnlineMask() const;
  [[nodiscard]] static bool areAllCpusOnline();
  void pinToFreeCPUs() const;
  /**
   * Pins the calling thread to @p mask. Unlike pinToFreeCPUs(), it calls no collective, so that a
   * thread of its own can call it with a mask that the main thread has computed.
   */
  static void pinCurrentThread(const async::as::CpuMask& mask);
  static std::string maskToString(const async::as::CpuMask& mask);
  static std::string maskToStringShort(const async::as::CpuMask& mask);
  [[nodiscard]] static async::as::CpuMask getNodeMask();
  void checkEnvVariables(utils::Env& env);
};

} // namespace seissol::parallel

#endif // SEISSOL_SRC_PARALLEL_PIN_H_
