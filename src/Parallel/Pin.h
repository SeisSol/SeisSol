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
#include <sched.h>
#include <string>

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
  /**
   * @brief The CPUs for threads that keep the halo exchange going on their own: those that its
   * libraries start (the proxies of NCCL and NVSHMEM), and the thread of actor-mpi.
   *
   * These are the CPUs of the OpenMP workers of this process without those of the calling thread:
   * the main thread polls all the time while the device computes, and the workers are idle then.
   * A thread that shares its CPU with a polling one gets it for a time slice at a time, which made
   * each message take milliseconds (Horizon, November 2026). If that leaves none, the free CPUs;
   * empty if there are none either.
   */
  [[nodiscard]] async::as::CpuMask getHelperMask() const;
  static bool freeCPUsMaskEmpty(const async::as::CpuMask& mask);
  /**
   * @return The lowest CPU in @p mask, or none if the mask is empty.
   */
  [[nodiscard]] static std::optional<int> firstCpu(const async::as::CpuMask& mask);
  [[nodiscard]] async::as::CpuMask getOnlineMask() const;
  [[nodiscard]] static bool areAllCpusOnline();
  void pinToFreeCPUs() const;
  static std::string maskToString(const async::as::CpuMask& mask);
  static std::string maskToStringShort(const async::as::CpuMask& mask);
  [[nodiscard]] static async::as::CpuMask getNodeMask();
  void checkEnvVariables();
};

/**
 * @brief Binds the calling thread to the CPUs of a mask while it exists, and back to its own ones
 * afterwards. The threads it starts meanwhile keep the mask: e.g. those of a library, which take
 * the CPUs of the thread that sets the library up. Does nothing for an empty mask.
 */
class ScopedAffinity {
  public:
  explicit ScopedAffinity(const async::as::CpuMask& mask);
  ~ScopedAffinity();

  ScopedAffinity(const ScopedAffinity&) = delete;
  ScopedAffinity& operator=(const ScopedAffinity&) = delete;
  ScopedAffinity(ScopedAffinity&&) = delete;
  ScopedAffinity& operator=(ScopedAffinity&&) = delete;

  private:
  bool bound_{false};
  async::as::CpuMask own_{};
};

} // namespace seissol::parallel

#endif // SEISSOL_SRC_PARALLEL_PIN_H_
