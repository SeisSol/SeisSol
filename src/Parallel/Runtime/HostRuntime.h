// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_HOSTRUNTIME_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_HOSTRUNTIME_H_

#include <algorithm>
#include <cstddef>
#include <functional>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace seissol::parallel::runtime {

/**
 * What the work of an action does with the data of its cluster, and with the data of the clusters
 * it observes (its neighbors). Besides, the work changes what only the cluster itself works on
 * (e.g. its DOFs).
 */
struct HostAccess {
  /// the work reads what the neighbors share with the cluster (e.g. their time data)
  bool readsNeighbors{false};

  /// the work writes what the cluster shares with its neighbors; it reads it in any case
  bool writesShared{false};
};

/**
 * The work of a cluster on the host: the counterpart of its StreamRuntime for the work of its
 * actions, and for the loops over its cells or faces therein.
 *
 * The work of an action runs right when the action is decided. A loop is a parallel region of its
 * own, which splits the range among all threads evenly, as a static schedule does.
 */
class HostRuntime {
  public:
  /**
   * Runs the work of an action. What it does with the data is stated by `access`.
   */
  template <typename F>
  void action(const HostAccess& /*access*/, F&& work) {
    std::invoke(std::forward<F>(work));
  }

  /**
   * Runs `body(begin, end)` on ranges that together cover [0, count) once each: one range per
   * thread, as a static schedule splits it. Whatever a range needs (e.g. scratch memory) can be set
   * up in the body, once per range.
   */
  template <typename F>
  static void forRange(std::size_t count, const F& body) {
    if (count == 0) {
      return;
    }
#ifdef _OPENMP
#pragma omp parallel default(none) shared(count, body)
    {
      const auto [begin, end] = staticRange(count,
                                            static_cast<std::size_t>(omp_get_thread_num()),
                                            static_cast<std::size_t>(omp_get_num_threads()));
      if (begin < end) {
        std::invoke(body, begin, end);
      }
    }
#else
    std::invoke(body, std::size_t{0}, count);
#endif
  }

  /**
   * Runs `body(i)` for every i in [0, count), split like forRange().
   */
  template <typename F>
  static void forEach(std::size_t count, const F& body) {
    forRange(count, [&body](std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {
        std::invoke(body, i);
      }
    });
  }

  /**
   * The range of thread `thread` of `threads` when [0, count) is split among them evenly: the
   * first `count % threads` threads take one element more than the others (as the static schedule
   * of LLVM's OpenMP runtime does).
   */
  static std::pair<std::size_t, std::size_t>
      staticRange(std::size_t count, std::size_t thread, std::size_t threads) {
    const auto base = count / threads;
    const auto extra = count % threads;
    const auto begin = thread * base + std::min(thread, extra);
    return {begin, begin + base + (thread < extra ? 1 : 0)};
  }
};

} // namespace seissol::parallel::runtime

#endif // SEISSOL_SRC_PARALLEL_RUNTIME_HOSTRUNTIME_H_
