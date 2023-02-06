#include "MiniSeisSolWriter.h"
#include "Parallel/MPI.h"
#include "utils/stringutils.h"
#include <filesystem>
#include <type_traits>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <unistd.h>

namespace seissol::writer::details {
template <typename T>
auto collect(T value) {
  auto collect = std::vector<T>(seissol::MPI::mpi.size());
  MPI_Datatype type;
  if constexpr (std::is_same_v<T, double>) {
    type = MPI_DOUBLE;
  } else if constexpr (std::is_same_v<T, int>) {
    type = MPI_INT;
  }

  MPI_Gather(&value, 1, type, collect.data(), 1, type, 0, seissol::MPI::mpi.comm());
  return collect;
}

auto collect(std::string& hostName) {
  auto localNameLength = static_cast<int>(hostName.size());
  auto lengths = collect(localNameLength);

  const size_t mpiSize = seissol::MPI::mpi.size();

  std::vector<int> displacements(mpiSize);
  std::exclusive_scan(lengths.begin(), lengths.end(), displacements.begin(), 0);

  const auto recvBufferSize = std::accumulate(lengths.begin(), lengths.end(), 0);
  std::vector<char> recvBuffer(recvBufferSize);

  MPI_Gatherv(hostName.c_str(),
              hostName.size(),
              MPI_CHAR,
              recvBuffer.data(),
              lengths.data(),
              displacements.data(),
              MPI_CHAR,
              0,
              seissol::MPI::mpi.comm());

  displacements.push_back(recvBufferSize);
  std::vector<std::string> hostNames(mpiSize);
  for (size_t rank = 0; rank < mpiSize; ++rank) {
    auto begin = displacements[rank];
    auto end = displacements[rank + 1];
    hostNames[rank] = std::string(&recvBuffer[begin], &recvBuffer[end]);
  }

  return hostNames;
}
} // namespace seissol::writer::details

void seissol::writer::MiniSeisSolWriter::write(double elapsedTime, double weight) {
  auto elapsedTimeVector = details::collect(elapsedTime);
  auto weightVector = details::collect(weight);

  auto localRanks = details::collect(seissol::MPI::mpi.sharedMemMpiRank());

  std::string hostName(256, ' ');
  if (gethostname(const_cast<char*>(hostName.c_str()), 1024) != 0) {
    hostName = "unknown-host";
  }
  utils::StringUtils::rtrim(hostName);
  auto hostNames = details::collect(hostName);

  if (seissol::MPI::mpi.rank() == 0) {
    std::vector<size_t> ranks(seissol::MPI::mpi.size());
    for (size_t i = 0; i < ranks.size(); ++i)
      ranks[i] = i;

    std::sort(ranks.begin(), ranks.end(), [&elapsedTimeVector](const size_t& i, const size_t& j) {
      return elapsedTimeVector[i] > elapsedTimeVector[j];
    });

    std::filesystem::path path(outputDirectory);
    path += std::filesystem::path("-miniSeissol.csv");

    std::fstream fileStream(path, std::ios::out);
    fileStream << "rank,localRank,hostname,elapsedTime,weight\n";

    for (auto rank : ranks) {
      fileStream << rank << ',' << localRanks[rank] << ',' << hostNames[rank] << ','
                 << elapsedTimeVector[rank] << ',' << weightVector[rank] << '\n';
    }

    fileStream.close();
  }
}