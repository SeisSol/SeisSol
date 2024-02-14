#include "MiniSeisSolWriter.h"
#include "Parallel/MPI.h"
#include "Common/filesystem.h"
#include <fstream>
#include <algorithm>
#include <unistd.h>

void seissol::writer::MiniSeisSolWriter::write(double elapsedTime, double weight) {
  auto elapsedTimeVector = seissol::MPI::mpi.collect(elapsedTime);
  auto weightVector = seissol::MPI::mpi.collect(weight);

  auto localRanks = seissol::MPI::mpi.collect(seissol::MPI::mpi.sharedMemMpiRank());

  if (seissol::MPI::mpi.rank() == 0) {
    std::vector<size_t> ranks(seissol::MPI::mpi.size());
    for (size_t i = 0; i < ranks.size(); ++i)
      ranks[i] = i;

    std::sort(ranks.begin(), ranks.end(), [&elapsedTimeVector](const size_t& i, const size_t& j) {
      return elapsedTimeVector[i] > elapsedTimeVector[j];
    });

    seissol::filesystem::path path(outputDirectory);
    path += seissol::filesystem::path("-miniSeissol.csv");

    std::fstream fileStream(path, std::ios::out);
    fileStream << "hostname,rank,localRank,elapsedTime,weight\n";

    const auto& hostNames = seissol::MPI::mpi.getHostNames();
    for (auto rank : ranks) {
      fileStream << "\"" << hostNames[rank] << "\"," << rank << ',' << localRanks[rank] << ','
                 << elapsedTimeVector[rank] << ',' << weightVector[rank] << '\n';
    }

    fileStream.close();
  }
}
