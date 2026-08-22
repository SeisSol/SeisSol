// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EasiReader.h"

#ifdef USE_EASI

#include "Reader/Datafield/AsagiReader.h"
#include "Reader/Scripting/DataTable.h"

#include <algorithm>
#include <cstddef>
#include <easi/Query.h>
#include <easi/ResultAdapter.h>
#include <easi/YAMLParser.h>
#include <easi/util/Slice.h>
#include <easi/util/Vector.h>
#include <exception>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utils/logger.h>
#include <vector>

namespace seissol::reader::scripting {

namespace {
// helper class to be independent from adapting single structs/arrays only
class MixedResultsAdapter : public easi::ResultAdapter {
  public:
  MixedResultsAdapter(std::size_t base,
                      const std::vector<DataEntry>& entries,
                      const std::vector<std::string>& provided)
      : base_(base), entries_(entries) {

    std::unordered_map<std::string, bool> providedMap;
    for (const auto& provide : provided) {
      providedMap[provide] = true;
    }

    for (std::size_t i = 0; i < entries.size(); ++i) {
      if (entries[i].direction != Direction::In) {
        indices_[entries[i].name] = i;
        if (providedMap.find(entries[i].name) == providedMap.end()) {
          logError() << "Easi script output parameter not found:" << entries[i].name;
        }
      }
    }
  }

  MixedResultsAdapter(std::size_t base,
                      const std::vector<DataEntry>& entries,
                      const std::vector<std::size_t>& indices)
      : base_(base), entries_(entries) {
    for (const auto& i : indices) {
      indices_[entries[i].name] = i;
    }
  }

  ~MixedResultsAdapter() override = default;
  void set(const std::string& parameter,
           const easi::Vector<unsigned>& index,
           const easi::Slice<double>& value) override {
    const auto& entryIdx = indices_.at(parameter);
    const auto& entry = entries_[entryIdx];

#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < value.size(); ++i) {
      const auto localIndex = base_ + index(i);
      entry.setValueAs(localIndex, value(i));
    }
  }
  [[nodiscard]] bool isSubset(const std::set<std::string>& parameters) const override {
    const auto myParams = this->parameters();
    return std::all_of(myParams.begin(), myParams.end(), [&](const auto& param) {
      return parameters.find(param) != parameters.end();
    });
  }
  ResultAdapter* subsetAdapter(const std::set<std::string>& subset) override {
    std::vector<std::size_t> subIndices;
    for (const auto& [_, index] : indices_) {
      if (subset.find(entries_[index].name) != subset.end()) {
        subIndices.emplace_back(index);
      }
    }
    return new MixedResultsAdapter(base_, entries_, subIndices);
  }
  [[nodiscard]] unsigned numberOfParameters() const override { return indices_.size(); }
  [[nodiscard]] std::set<std::string> parameters() const override {
    std::set<std::string> params;
    for (const auto& [_, index] : indices_) {
      params.insert(entries_[index].name);
    }
    return params;
  }

  private:
  std::size_t base_{};
  const std::vector<DataEntry>& entries_;
  std::unordered_map<std::string, std::size_t> indices_;
};
} // namespace

EasiReader::~EasiReader() = default;

EasiReader::EasiReader(const std::string& script, const std::vector<std::string>& inVars)
    : script_(script) {
#ifdef USE_ASAGI
  asagiReader_ = std::make_unique<seissol::asagi::AsagiReader>();
#else
  asagiReader_.reset();
#endif

  const auto dimensionNames = std::set<std::string>(inVars.begin(), inVars.end());
  inVars_ = std::vector<std::string>(dimensionNames.begin(), dimensionNames.end());
  parser_ = std::make_unique<easi::YAMLParser>(dimensionNames, asagiReader_.get());

  try {
    components_ = std::unique_ptr<easi::Component>(parser_->parse(script));
  } catch (const std::exception& error) {
    logError() << "Error while parsing easi file" << script << ":" << std::string(error.what());
  }

  const auto outVarsPre = components_->suppliedParameters();
  outVars_ = std::vector<std::string>(outVarsPre.begin(), outVarsPre.end());
}

void EasiReader::call(const scripting::DataTable& table) {
  // avoid overflows for old easi versions
  constexpr std::size_t CallBatchSize = 1ULL << 31;

  const std::size_t batches = (table.numPoints() + CallBatchSize - 1) / CallBatchSize;

  const auto& entries = table.dataEntries();

  // a "magic" entry for easi; will _not_ be considered as input/output
  std::optional<std::size_t> groupEntry;

  std::unordered_map<std::string, std::size_t> inVarMap;
  for (std::size_t i = 0; i < inVars_.size(); ++i) {
    inVarMap[inVars_[i]] = i;
  }
  std::vector<std::size_t> inVarToEntry(inVars_.size());
  std::vector<bool> found(inVars_.size());

  for (std::size_t i = 0; i < entries.size(); ++i) {
    const auto& entry = entries[i];
    if (entry.name == "group") {
      // special behavior: group is int
      groupEntry = i;
    }
    if (inVarMap.find(entry.name) != inVarMap.end()) {
      const auto index = inVarMap.at(entry.name);
      inVarToEntry[index] = i;
      found[index] = true;
    }
  }

  for (std::size_t i = 0; i < inVars_.size(); ++i) {
    if (!found[i]) {
      logError() << "Error while loading easi script input parameters:" << inVars_[i]
                 << "not found.";
    }
  }

  easi::Query query(0, 0);
  for (std::size_t batch = 0; batch < batches; ++batch) {
    const auto base = CallBatchSize * batch;
    const auto size = std::min(CallBatchSize, table.numPoints() - base);
    if (query.index.size() != size) {
      // reallocate
      query = easi::Query(size, inVars_.size());
    }

    auto adapter = MixedResultsAdapter(batch * CallBatchSize, entries, outVars_);

#pragma omp parallel for schedule(static)
    for (std::size_t point = 0; point < size; ++point) {
      for (std::size_t i = 0; i < inVars_.size(); ++i) {
        const auto& entry = entries[inVarToEntry[i]];
        query.x(point, i) = entry.getValueAs<double>(base + point);
      }
      if (groupEntry.has_value()) {
        const auto& entry = entries[groupEntry.value()];
        query.group(point) = entry.getValueAs<int32_t>(base + point);
      }
    }

    try {
      components_->evaluate(query, adapter);
    } catch (const std::exception& error) {
      logError() << "Error while applying easi file" << script_ << ":" << std::string(error.what());
    }
  }
}

} // namespace seissol::reader::scripting

#endif
