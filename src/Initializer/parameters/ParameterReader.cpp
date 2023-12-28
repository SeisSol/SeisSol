#include "ParameterReader.h"

namespace seissol::initializers::parameters {

void ParameterReader::warnDeprecatedSingle(const std::string& field) {
  if (hasField(field)) {
    visited.emplace(field);
    logInfo(seissol::MPI::mpi.rank())
        << "The field" << field
        << "is no longer in use. You may safely remove it from your parameters file.";
  }
}

void ParameterReader::warnDeprecated(const std::vector<std::string>& fields) {
  for (const auto& field : fields) {
    warnDeprecatedSingle(field);
  }
}

void ParameterReader::warnUnknown(const std::string& prefix) const {
  for (const auto& subnodes : node) {
    auto field = subnodes.first.as<std::string>();
    if (visited.find(field) == visited.end()) {
      logWarning(seissol::MPI::mpi.rank())
          << "The field" << field << "in" << prefix
          << "was given in the parameter file, but is unknown to SeisSol.";
    }
  }
  for (const auto& pair : subreaders) {
    pair.second.warnUnknown(pair.first);
  }
}

void ParameterReader::markUnused(const std::vector<std::string>& fields) {
  for (const auto& field : fields) {
    logDebug(seissol::MPI::mpi.rank()) << "The field" << field << "is ignored (if it is found).";
    visited.emplace(field);
  }
}

ParameterReader& ParameterReader::readSubNode(const std::string& subnodeName) {
  visited.emplace(subnodeName);
  logDebug(seissol::MPI::mpi.rank()) << "Entering section" << subnodeName;
  if (subreaders.find(subnodeName) == subreaders.end()) {
    if (hasField(subnodeName)) {
      subreaders.emplace(subnodeName, ParameterReader(node[subnodeName], false));
    } else {
      logDebug(seissol::MPI::mpi.rank())
          << "Section" << subnodeName
          << "not found in the given parameter file. Using an empty reader.";
      subreaders.emplace(subnodeName, ParameterReader(node[subnodeName], true));
    }
  }
  return subreaders.at(subnodeName);
}

bool ParameterReader::hasField(const std::string& field) { return !empty && node[field]; }

} // namespace seissol::initializers::parameters
