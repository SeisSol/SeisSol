// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MeasuredWeights.h"

#include "Geometry/MeshReader.h"
#include "IO/Reader/Distribution.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Clustering/VertexWeights/WeightsModels.h"
#include "Initializer/MemoryManager.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Colormap.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/LoopStatistics.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <hdf5.h>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

namespace {

static_assert(sizeof(std::size_t) == sizeof(std::uint64_t), "the cell ids go to the file as such");

//! stops with a message, if an HDF5 call failed
template <typename T>
T checked(T status, const char* what, const std::string& file) {
  if (status < 0) {
    logError() << "Writing the measured cell costs to" << file << "failed in" << what;
  }
  return status;
}

/**
 * Where the cells of this rank start in a sequence of the cells of all ranks, one after the other;
 * and how many there are in all. Collective.
 */
std::pair<std::uint64_t, std::uint64_t> runOf(std::uint64_t count, MPI_Comm comm) {
  std::uint64_t first = 0;
  MPI_Exscan(&count, &first, 1, MPI_UINT64_T, MPI_SUM, comm);
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  // MPI_Exscan leaves the result on the first rank undefined
  if (rank == 0) {
    first = 0;
  }
  std::uint64_t total = 0;
  MPI_Allreduce(&count, &total, 1, MPI_UINT64_T, MPI_SUM, comm);
  return {first, total};
}

//! writes `data` of every rank to the dataset `name`, one rank after the other
void writeRun(hid_t file,
              const char* name,
              hid_t type,
              const void* data,
              std::uint64_t first,
              std::uint64_t count,
              std::uint64_t total,
              const std::string& fileName) {
  const hsize_t globalSize = total;
  const hid_t fileSpace = checked(H5Screate_simple(1, &globalSize, nullptr), name, fileName);
  const hsize_t localSize = count;
  const hid_t memorySpace = checked(H5Screate_simple(1, &localSize, nullptr), name, fileName);
  if (count > 0) {
    const hsize_t start = first;
    checked(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, &start, nullptr, &localSize, nullptr),
            name,
            fileName);
  } else {
    checked(H5Sselect_none(fileSpace), name, fileName);
    checked(H5Sselect_none(memorySpace), name, fileName);
  }
  const hid_t dataset =
      checked(H5Dcreate(file, name, type, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
              name,
              fileName);
  const hid_t transfer = checked(H5Pcreate(H5P_DATASET_XFER), name, fileName);
  checked(H5Pset_dxpl_mpio(transfer, H5FD_MPIO_COLLECTIVE), name, fileName);
  checked(H5Dwrite(dataset, type, memorySpace, fileSpace, transfer, data), name, fileName);
  checked(H5Pclose(transfer), name, fileName);
  checked(H5Dclose(dataset), name, fileName);
  checked(H5Sclose(memorySpace), name, fileName);
  checked(H5Sclose(fileSpace), name, fileName);
}

} // namespace

std::optional<MeasuredCellCosts> measureCellCosts(const LoopStatistics& statistics,
                                                  MemoryManager& memoryManager,
                                                  const geometry::MeshReader& meshReader,
                                                  double simulatedTime,
                                                  MPI_Comm comm) {
  const auto cellRegions = std::array{statistics.getRegion("computeLocalIntegration"),
                                      statistics.getRegion("computeNeighboringIntegration"),
                                      statistics.getRegion("computePointSources")};
  const auto faceRegion = statistics.getRegion("computeDynamicRupture");

  auto& ltsStorage = memoryManager.ltsStorage();
  auto& drStorage = memoryManager.drStorage();
  const auto& colorMap = ltsStorage.getColorMap();
  const auto colors = colorMap.size();

  // The time of the cells goes by group: a time cluster in a configuration, its interior and copy
  // layer together. Per layer, the cells at the boundary of the partition would cost more than
  // the others (the copy layers are the smaller ones, so their fixed costs per loop count for
  // fewer cells; and some copy cells get computed more than once), which says nothing about where
  // the boundary of the next partition should go. The layers of faces are laid out along the same
  // colors as those of the cells.
  const auto groupOf = [&](std::size_t color) {
    auto identifier = ltsStorage.layer(color).getIdentifier();
    identifier.halo = HaloType::Interior;
    return colorMap.colorId(identifier);
  };
  std::vector<double> cellTime(colors, 0.0);
  std::vector<double> faceTime(colors, 0.0);
  double measured = 0;
  for (std::size_t color = 0; color < colors; ++color) {
    if (ltsStorage.layer(color).getIdentifier().halo == HaloType::Ghost) {
      continue;
    }
    double time = 0;
    for (const auto region : cellRegions) {
      time += statistics.layerTime(region, color);
    }
    const auto face = statistics.layerTime(faceRegion, color);
    cellTime[groupOf(color)] += time;
    faceTime[groupOf(color)] += face;
    measured += time + face;
  }

  int timed = measured > 0 ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &timed, 1, MPI_INT, MPI_MIN, comm);
  if (timed == 0) {
    return {};
  }

  const auto& elements = meshReader.getElements();
  const auto& backmap = memoryManager.backmap();
  const auto& faults = meshReader.getFault();
  const auto& drBackmap = memoryManager.drBackmap();

  // the cells and faces of each group, each once
  std::vector<std::size_t> cellGroup(elements.size());
  std::vector<double> cells(colors, 0.0);
  for (std::size_t cell = 0; cell < elements.size(); ++cell) {
    cellGroup[cell] = groupOf(backmap.get(cell).color);
    ++cells[cellGroup[cell]];
  }
  std::vector<double> faces(colors, 0.0);
  for (std::size_t face = 0; face < faults.size(); ++face) {
    const auto position = drBackmap.getDup(face, 0);
    if (position.has_value()) {
      ++faces[groupOf(position->color)];
    }
  }

  MeasuredCellCosts costs;
  costs.ids.resize(elements.size());
  costs.costs.resize(elements.size());
  for (std::size_t cell = 0; cell < elements.size(); ++cell) {
    costs.ids[cell] = elements[cell].globalId;
    costs.costs[cell] = cellTime[cellGroup[cell]] / cells[cellGroup[cell]];
  }

  // a face goes to the cells at it that are on this rank; a face between two ranks is computed on
  // both
  for (std::size_t face = 0; face < faults.size(); ++face) {
    const auto position = drBackmap.getDup(face, 0);
    if (!position.has_value()) {
      continue;
    }
    const auto group = groupOf(position->color);
    const auto& fault = faults[face];
    const bool plus = fault.element.hasValue();
    const bool minus = fault.neighborElement.hasValue();
    const double share = plus && minus ? 0.5 : 1.0;
    const double time = share * faceTime[group] / faces[group];
    if (plus) {
      costs.costs[fault.element.value()] += time;
    }
    if (minus) {
      costs.costs[fault.neighborElement.value()] += time;
    }
  }

  if (simulatedTime > 0) {
    for (auto& cost : costs.costs) {
      cost /= simulatedTime;
    }
  }
  return costs;
}

void writeCellCosts(const std::string& file,
                    const MeasuredCellCosts& costs,
                    double simulatedTime,
                    MPI_Comm comm) {
  const auto [first, total] = runOf(costs.ids.size(), comm);

  const hid_t access = checked(H5Pcreate(H5P_FILE_ACCESS), "H5Pcreate", file);
  checked(H5Pset_fapl_mpio(access, comm, MPI_INFO_NULL), "H5Pset_fapl_mpio", file);
  const hid_t h5file =
      checked(H5Fcreate(file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access), "H5Fcreate", file);
  checked(H5Pclose(access), "H5Pclose", file);

  writeRun(
      h5file, "ids", H5T_NATIVE_UINT64, costs.ids.data(), first, costs.ids.size(), total, file);
  writeRun(h5file,
           "costs",
           H5T_NATIVE_DOUBLE,
           costs.costs.data(),
           first,
           costs.costs.size(),
           total,
           file);

  // the simulated time the costs were measured over, for whoever looks at them
  const hid_t scalar = checked(H5Screate(H5S_SCALAR), "H5Screate", file);
  const hid_t attribute = checked(
      H5Acreate(h5file, "simulatedTime", H5T_NATIVE_DOUBLE, scalar, H5P_DEFAULT, H5P_DEFAULT),
      "H5Acreate",
      file);
  checked(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &simulatedTime), "H5Awrite", file);
  checked(H5Aclose(attribute), "H5Aclose", file);
  checked(H5Sclose(scalar), "H5Sclose", file);
  checked(H5Fclose(h5file), "H5Fclose", file);
}

std::vector<double> distributeCellCosts(const MeasuredCellCosts& costs,
                                        std::size_t firstCell,
                                        std::size_t count,
                                        MPI_Comm comm) {
  const auto [measuredFirst, measuredTotal] = runOf(costs.ids.size(), comm);
  const auto [first, total] = runOf(count, comm);
  if (measuredTotal != total) {
    logError() << "The measured costs are those of" << measuredTotal << "cells, but the mesh has"
               << total << "cells.";
  }

  std::vector<std::size_t> cells(count);
  std::iota(cells.begin(), cells.end(), firstCell);
  io::reader::Distributor distributor(comm);
  distributor.setup(costs.ids, cells);
  std::vector<double> result(count);
  distributor.distribute(result.data(), costs.costs.data()).complete();
  return result;
}

std::vector<double> readCellCosts(const std::string& file,
                                  std::size_t firstCell,
                                  std::size_t count,
                                  MPI_Comm comm) {
  MeasuredCellCosts costs;
  io::reader::file::Hdf5Reader reader(comm);
  reader.openFile(file);
  costs.ids = reader.readData<std::size_t>("ids");
  costs.costs = reader.readData<double>("costs");
  reader.closeFile();
  if (costs.ids.size() != costs.costs.size()) {
    logError() << "The file" << file << "holds" << costs.ids.size() << "cell ids, but"
               << costs.costs.size() << "costs.";
  }
  return distributeCellCosts(costs, firstCell, count, comm);
}

MeasuredWeights::MeasuredWeights(std::string file) : file_(std::move(file)) {}

MeasuredWeights::MeasuredWeights(const MeasuredCellCosts& costs) : costs_(&costs) {}

void MeasuredWeights::setVertexWeights() {
  const auto& clustering = this->clustering();
  const auto comm = seissol::Mpi::mpi.comm();
  const std::size_t count = clustering.clusterIds.size();

  // the cells of this rank are a run of those of the file
  const auto [first, total] = runOf(count, comm);
  if (costs_ != nullptr) {
    logInfo() << "Partition by the costs of the cells measured in the run before.";
  } else {
    logInfo() << "Partition by the costs of the cells measured in" << file_;
  }
  const auto costs = costs_ != nullptr ? distributeCellCosts(*costs_, first, count, comm)
                                       : readCellCosts(file_, first, count, comm);

  // the weights sum to what those of the default model would, which ParMETIS has been taking
  ExponentialWeights reference;
  reference.build(clustering);
  std::array<double, 2> sums{0, 0};
  for (std::size_t cell = 0; cell < count; ++cell) {
    sums[0] += static_cast<double>(reference.vertexWeights()[cell]);
    sums[1] += costs[cell];
  }
  MPI_Allreduce(MPI_IN_PLACE, sums.data(), sums.size(), MPI_DOUBLE, MPI_SUM, comm);
  if (!(sums[1] > 0)) {
    logError() << "The measured costs of the cells sum to" << sums[1];
  }
  const double scale = sums[0] / sums[1];

  std::array<double, 2> range{std::numeric_limits<double>::max(), 0};
  for (std::size_t cell = 0; cell < count; ++cell) {
    const auto weight = std::max(1.0, std::round(costs[cell] * scale));
    vertexWeights_[cell] = static_cast<std::uint64_t>(weight);
    range[0] = std::min(range[0], weight);
    range[1] = std::max(range[1], weight);
  }
  range[0] = -range[0];
  MPI_Allreduce(MPI_IN_PLACE, range.data(), range.size(), MPI_DOUBLE, MPI_MAX, comm);
  const double mean = sums[0] / static_cast<double>(total);
  logInfo() << "The weights of the cells range from" << -range[0] / mean << "to" << range[1] / mean
            << "times their mean.";
}

void MeasuredWeights::setAllowedImbalances() {
  // as for ExponentialWeights
  constexpr double TinyLtsWeightImbalance{1.01};
  imbalances_.assign(1, TinyLtsWeightImbalance);
}

} // namespace seissol::initializer
