#pragma once

#include <algorithm>
#include <optional>
#include <unordered_map>
#include <vector>

namespace seissol::initializer {

struct ClusterLayout {
int rate;
double minimumTimestep;
std::size_t globalClusterCount;
std::vector<int> localClusterIds;
std::unordered_map<int, std::size_t> localClusterIdsInverse;

ClusterLayout(int rate, double minimumTimestep, std::size_t globalClusterCount, const std::vector<int>& localClusterIds) : rate(rate), minimumTimestep(minimumTimestep), globalClusterCount(globalClusterCount), localClusterIds(localClusterIds) {
    std::sort(this->localClusterIds.begin(), this->localClusterIds.end());
    for (std::size_t i = 0; i < this->localClusterIds.size(); ++i) {
        localClusterIdsInverse[this->localClusterIds[i]] = i;
    }
}

std::optional<std::size_t> localClusterFromGlobal(int id) {
    if (localClusterIdsInverse.find(id) != localClusterIdsInverse.end()) {
        return id;
    }
    else {
        return {};
    }
}

double timestepRate(std::size_t id) {
    return clusterRate(id) * minimumTimestep;
}

long clusterRate(std::size_t id) {
    long value = 1;
    const int cluster = localClusterIds.at(id);
    for (int i = 0; i < cluster; ++i) {
        value *= rate;
    }
    return value;
}

};

} // namespace seissol::initializer
