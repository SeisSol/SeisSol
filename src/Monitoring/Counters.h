// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_COUNTERS_H_
#define SEISSOL_SRC_MONITORING_COUNTERS_H_

#include "Monitoring/Unit.h"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace seissol::monitoring {

//! @brief How the values that a counter has on the ranks make up its value for the run.
enum class CounterScope {
  //! Every rank counts what it did, and the run did the sum over the ranks.
  Rank,
  //! One rank of a node counts for the whole node and the others count nothing, so that the run
  //! is the sum over the nodes.
  Node
};

//! @brief What a counter counts, and what the summary and the tables call it and its rate.
struct CounterInfo {
  //! e.g. hardwareFlop or nodeEnergy
  std::string name;
  //! e.g. hardwareFlopPerSecond or nodePower
  std::string rateName;
  const SIUnit* unit{&UnitScalar};
  const SIUnit* rateUnit{&UnitScalar};
  CounterScope scope{CounterScope::Rank};
  //! Whether the summary gives it for each stage of the run, which is worth it for what the
  //! initialization takes as well, such as energy, and not for what only the simulation does.
  bool perStage{false};
};

/**
 * @brief A value of a counter: what it counted since it was registered.
 *
 * A counter that keeps a clock of its own, such as an energy sensor, gives the time of the value on
 * that clock, in seconds, which its rate is taken over; the others take the time of the rank.
 */
struct CounterValue {
  double amount{0};
  std::optional<double> time;
};

/**
 * @brief The counters of a run: quantities that only grow, such as FLOP, cell updates or energy.
 *
 * Every rank registers the same counters, in the same order. A rank that has nothing to count for a
 * counter gives no value for it, which is how a counter of a node is counted once per node. A
 * snapshot reads all counters on a rank; what they counted over a stretch of the run is the
 * difference of two snapshots, summed over the ranks.
 */
class Counters {
  public:
  using Reader = std::function<std::optional<CounterValue>()>;

  /**
   * @brief Registers a counter, which @p reader reads; returns its position.
   *
   * A counter counts from zero when it is registered: a snapshot taken before counts as zero for
   * it.
   */
  std::size_t add(CounterInfo info, Reader reader);

  [[nodiscard]] std::size_t size() const { return counters_.size(); }

  [[nodiscard]] const CounterInfo& info(std::size_t counter) const;

  //! @brief The position of the counter named @p name, if there is one.
  [[nodiscard]] std::optional<std::size_t> find(std::string_view name) const;

  //! @brief The values of the counters on this rank, at one moment.
  struct Snapshot {
    //! the time of the rank, in seconds
    double time{0};
    //! one per counter registered at the moment; none where the rank gave no value
    std::vector<std::optional<CounterValue>> values;
  };

  //! @brief Reads all counters on this rank; @p time is the time of the rank, in seconds.
  [[nodiscard]] Snapshot snapshot(double time) const;

  //! @brief What the counters counted over a stretch of the run, over all ranks.
  struct Totals {
    //! the difference of each counter, summed over the ranks
    std::vector<double> amounts;
    //! the difference of each counter over the time of the rank (or of its own clock), summed over
    //! the ranks
    std::vector<double> rates;
    //! the number of ranks that gave a value for the counter at both ends of the stretch: the
    //! nodes for a counter of the nodes
    std::vector<int> contributors;

    //! @brief Whether the counter was counted by any rank at all.
    [[nodiscard]] bool counted(std::size_t counter) const {
      return counter < contributors.size() && contributors[counter] > 0;
    }
  };

  /**
   * @brief What the counters counted between @p from and @p to, summed over all ranks. Collective.
   *
   * The result is the same on all ranks.
   */
  [[nodiscard]] Totals reduce(const Snapshot& from, const Snapshot& to) const;

  /**
   * @brief The same as reduce(), for several stretches at once, with a single reduction.
   * Collective.
   */
  [[nodiscard]] std::vector<Totals>
      reduce(const std::vector<std::pair<const Snapshot*, const Snapshot*>>& stretches) const;

  private:
  struct Counter {
    CounterInfo info;
    Reader reader;
  };

  std::vector<Counter> counters_;
};

} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_COUNTERS_H_
