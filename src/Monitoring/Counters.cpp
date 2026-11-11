// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Counters.h"

#include "Parallel/MPI.h"

#include <cmath>
#include <cstddef>
#include <mpi.h>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

namespace seissol::monitoring {

std::size_t Counters::add(CounterInfo info, Reader reader) {
  counters_.push_back(Counter{.info = std::move(info), .reader = std::move(reader)});
  return counters_.size() - 1;
}

const CounterInfo& Counters::info(std::size_t counter) const { return counters_.at(counter).info; }

std::optional<std::size_t> Counters::find(std::string_view name) const {
  for (std::size_t i = 0; i < counters_.size(); ++i) {
    if (counters_[i].info.name == name) {
      return i;
    }
  }
  return {};
}

Counters::Snapshot Counters::snapshot(double time) const {
  Snapshot snapshot;
  snapshot.time = time;
  snapshot.values.reserve(counters_.size());
  for (const auto& counter : counters_) {
    snapshot.values.push_back(counter.reader());
  }
  return snapshot;
}

Counters::Totals Counters::reduce(const Snapshot& from, const Snapshot& to) const {
  return reduce({{&from, &to}}).front();
}

std::vector<Counters::Totals> Counters::reduce(
    const std::vector<std::pair<const Snapshot*, const Snapshot*>>& stretches) const {
  // per stretch and counter: the difference, the rate and whether the rank counted it
  constexpr std::size_t Entries = 3;
  const auto count = counters_.size();
  std::vector<double> buffer(stretches.size() * count * Entries, 0.0);

  // a counter registered after a snapshot counted nothing before it
  const auto value = [](const Snapshot& snapshot,
                        std::size_t counter) -> std::optional<CounterValue> {
    if (counter >= snapshot.values.size()) {
      return CounterValue{};
    }
    return snapshot.values[counter];
  };

  for (std::size_t stretch = 0; stretch < stretches.size(); ++stretch) {
    const auto& [from, to] = stretches[stretch];
    for (std::size_t counter = 0; counter < count; ++counter) {
      const auto start = value(*from, counter);
      const auto end = value(*to, counter);
      if (!start.has_value() || !end.has_value()) {
        continue;
      }
      const double amount = end->amount - start->amount;
      const double time = start->time.has_value() && end->time.has_value()
                              ? end->time.value() - start->time.value()
                              : to->time - from->time;
      auto* entry = &buffer[(stretch * count + counter) * Entries];
      entry[0] = amount;
      entry[1] = time > 0 ? amount / time : 0.0;
      entry[2] = 1.0;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,
                buffer.data(),
                static_cast<int>(buffer.size()),
                MPI_DOUBLE,
                MPI_SUM,
                seissol::Mpi::mpi.comm());

  std::vector<Totals> totals(stretches.size());
  for (std::size_t stretch = 0; stretch < stretches.size(); ++stretch) {
    auto& result = totals[stretch];
    result.amounts.resize(count);
    result.rates.resize(count);
    result.contributors.resize(count);
    for (std::size_t counter = 0; counter < count; ++counter) {
      const auto* entry = &buffer[(stretch * count + counter) * Entries];
      result.amounts[counter] = entry[0];
      result.rates[counter] = entry[1];
      result.contributors[counter] = static_cast<int>(std::lround(entry[2]));
    }
  }
  return totals;
}

} // namespace seissol::monitoring
