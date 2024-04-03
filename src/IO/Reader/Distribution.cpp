#include "Distribution.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <unordered_map>
#include <vector>

namespace {
static int getRank(std::size_t id, std::size_t count, int commsize) {
  const std::size_t piece1 = (count + commsize - 1) / commsize;
  const std::size_t rest = count % commsize;
  const std::size_t part1 = piece1 * rest;
  if (id < part1) {
    return id / piece1;
  } else {
    const std::size_t piece2 = count / commsize;
    return (id - part1) / piece2 + rest;
  }
}

template <typename T>
static std::vector<T> computeHistogram(const std::vector<T>& input) {
  std::vector<T> output(input.size() + 1);
  for (std::size_t i = 1; i < output.size(); ++i) {
    output[i] = output[i - 1] + input[i - 1];
  }
  return output;
}

template <typename T>
static std::vector<std::pair<T, int>> distributeIds(const std::vector<std::pair<T, int>>& source,
                                                    MPI_Comm comm,
                                                    MPI_Datatype sizetype,
                                                    MPI_Datatype datatype,
                                                    int tag) {
  int commsize;
  MPI_Comm_size(comm, &commsize);

  std::vector<std::size_t> sourceToSend(commsize);
  for (const auto& [_, j] : source) {
    ++sourceToSend[j];
  }
  std::vector<std::size_t> sourceToRecv(commsize);

  MPI_Alltoall(sourceToSend.data(), 1, sizetype, sourceToRecv.data(), 1, sizetype, comm);

  auto sourceToSendHistogram = computeHistogram(sourceToSend);

  std::vector<T> sortedSource(source.size());
  {
    auto tempHistogram = sourceToSendHistogram;
    for (const auto& [i, j] : source) {
      auto& position = tempHistogram[j];
      sortedSource[position] = i;
      ++position;
    }
  }

  auto sourceToRecvHistogram = computeHistogram(sourceToRecv);

  std::vector<T> intermediate(sourceToRecvHistogram.back());

  {
    std::vector<MPI_Request> requests(static_cast<std::size_t>(commsize) * 2, MPI_REQUEST_NULL);

    for (int i = 0; i < commsize; ++i) {
      if (sourceToSend[i] > 0) {
        MPI_Isend(sortedSource.data() + sourceToSendHistogram[i],
                  sourceToSend[i],
                  datatype,
                  i,
                  tag,
                  comm,
                  &requests[i]);
      }
      if (sourceToRecv[i] > 0) {
        MPI_Irecv(intermediate.data() + sourceToRecvHistogram[i],
                  sourceToRecv[i],
                  datatype,
                  i,
                  tag,
                  comm,
                  &requests[static_cast<std::size_t>(commsize) + i]);
      }
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  }

  std::vector<std::pair<T, int>> outputIntermediate(intermediate.size());
  int rank = 0;
  for (std::size_t i = 0; i < intermediate.size(); ++i) {
    if (sourceToRecvHistogram[rank + 1] <= i) {
      ++rank;
    }
    outputIntermediate[i] = std::pair<T, int>(intermediate[i], rank);
  }

  return outputIntermediate;
}

static std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
    matchRanks(const std::vector<std::pair<std::size_t, int>>& intermediateSource,
               const std::vector<std::pair<std::size_t, int>>& intermediateTarget,
               const std::vector<std::size_t>& source,
               MPI_Comm comm,
               MPI_Datatype sizetype,
               MPI_Datatype datatype,
               int tag) {
  int commsize;
  MPI_Comm_size(comm, &commsize);

  std::vector<std::size_t> sendOffsets, sendReorder;

  std::vector<std::pair<std::pair<std::size_t, int>, int>> sourceToTargetRankMap(
      intermediateSource.size());
  {
    std::unordered_multimap<std::size_t, int> intermediateSourceMap;
    for (const auto& i : intermediateSource) {
      intermediateSourceMap.insert(i);
    }
    for (std::size_t i = 0; i < intermediateSource.size(); ++i) {
      for (auto it = intermediateSourceMap.find(intermediateTarget[i].first);
           it != intermediateSourceMap.end() && it->first == intermediateTarget[i].first;
           ++it) { // TODO: handle the multimap properly
        sourceToTargetRankMap[i] = {{it->first, intermediateTarget[i].second}, it->second};
      }
    }
  }

  auto sourceToTarget = distributeIds<std::pair<std::size_t, int>>(
      sourceToTargetRankMap, comm, sizetype, datatype, tag);
  {
    std::unordered_map<std::size_t, std::size_t> sourceMap;
    for (std::size_t i = 0; i < source.size(); ++i) {
      sourceMap[source[i]] = i;
    }

    // TODO: we already have that: sourceToRecvHistogram
    std::vector<std::size_t> sourceCounter(commsize);
    for (std::size_t i = 0; i < source.size(); ++i) {
      ++sourceCounter[sourceToTarget[i].first.second];
    }
    sendOffsets = computeHistogram(sourceCounter);
    sendReorder.resize(sourceToTarget.size());
    for (std::size_t i = 0; i < sourceToTarget.size(); ++i) {
      sendReorder[i] = sourceMap.at(sourceToTarget[i].first.first);
    }
  }

  return {sendOffsets, sendReorder};
}
} // namespace

namespace seissol::io::reader {
Distributor::Distributor(MPI_Comm comm) : comm(comm) {}

void Distributor::setup(const std::vector<std::size_t>& sourceIds,
                        const std::vector<std::size_t>& targetIds) {
  int commsize;
  int commrank;

  constexpr int TagToIntermediateSource = 20;
  constexpr int TagToIntermediateTarget = 21;
  constexpr int TagFromIntermediateSource = 22;
  constexpr int TagFromIntermediateTarget = 23;

  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &commrank);

  MPI_Datatype sizetype =
      seissol::io::datatype::convertToMPI(seissol::io::datatype::inferDatatype<std::size_t>());
  MPI_Datatype pairtype = seissol::io::datatype::convertToMPI(
      seissol::io::datatype::inferDatatype<std::pair<std::size_t, int>>());

  auto sourceMax = *std::max_element(sourceIds.begin(), sourceIds.end());
  auto targetMax = *std::max_element(targetIds.begin(), targetIds.end());
  std::size_t localMax = std::max(sourceMax, targetMax);
  std::size_t globalMax;
  MPI_Allreduce(&localMax, &globalMax, 1, sizetype, MPI_MAX, comm);
  auto globalCount = globalMax + 1;

  // note that the following operations all need to be stable (i.e. order-preserving) for these
  // methods to work correctly right now

  std::vector<std::pair<std::size_t, int>> source(sourceIds.size());
  for (std::size_t i = 0; i < source.size(); ++i) {
    source[i] =
        std::pair<std::size_t, int>(sourceIds[i], getRank(sourceIds[i], globalCount, commsize));
  }
  std::vector<std::pair<std::size_t, int>> target(targetIds.size());
  for (std::size_t i = 0; i < target.size(); ++i) {
    target[i] =
        std::pair<std::size_t, int>(targetIds[i], getRank(targetIds[i], globalCount, commsize));
  }

  auto intermediateSource =
      distributeIds<std::size_t>(source, comm, sizetype, sizetype, TagToIntermediateSource);
  auto intermediateTarget =
      distributeIds<std::size_t>(target, comm, sizetype, sizetype, TagToIntermediateTarget);

  auto sendResult = matchRanks(intermediateSource,
                               intermediateTarget,
                               sourceIds,
                               comm,
                               sizetype,
                               pairtype,
                               TagFromIntermediateSource);
  sendOffsets = sendResult.first;
  sendReorder = sendResult.second;

  auto recvResult = matchRanks(intermediateSource,
                               intermediateTarget,
                               targetIds,
                               comm,
                               sizetype,
                               pairtype,
                               TagFromIntermediateTarget);
  recvOffsets = recvResult.first;
  recvReorder = recvResult.second;
}

void Distributor::distributeInternal(void* target, const void* source, MPI_Datatype datatype) {
  constexpr int Tag = 30;

  int commsize;
  int commrank;

  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &commrank);

  int typesizeInt;
  MPI_Type_size(datatype, &typesizeInt);
  std::size_t typesize = typesizeInt;

  std::vector<MPI_Request> requests(static_cast<std::size_t>(commsize) * 2, MPI_REQUEST_NULL);

  const char* sourceChar = reinterpret_cast<const char*>(source);
  char* targetChar = reinterpret_cast<char*>(target);

  char* sourceReordered = reinterpret_cast<char*>(std::malloc(typesize * sendReorder.size()));
  char* targetReordered = reinterpret_cast<char*>(std::malloc(typesize * recvReorder.size()));

  for (std::size_t i = 0; i < sendReorder.size(); ++i) {
    std::memcpy(sourceReordered + i * typesize, sourceChar + sendReorder[i] * typesize, typesize);
  }

  for (int i = 0; i < commsize; ++i) {
    if (sendOffsets[i + 1] > sendOffsets[i]) {
      MPI_Isend(sourceReordered + sendOffsets[i] * typesize,
                sendOffsets[i + 1] - sendOffsets[i],
                datatype,
                i,
                Tag,
                comm,
                &requests[i]);
    }
    if (recvOffsets[i + 1] > recvOffsets[i]) {
      MPI_Irecv(targetReordered + recvOffsets[i] * typesize,
                recvOffsets[i + 1] - recvOffsets[i],
                datatype,
                i,
                Tag,
                comm,
                &requests[static_cast<std::size_t>(commsize) + i]);
    }
  }
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

  for (std::size_t i = 0; i < recvReorder.size(); ++i) {
    std::memcpy(targetChar + recvReorder[i] * typesize, targetReordered + i * typesize, typesize);
  }

  std::free(sourceReordered);
  std::free(targetReordered);
}
} // namespace seissol::io::reader
