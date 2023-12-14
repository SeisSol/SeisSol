// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "Args.hpp"
#include "LtsStructureParser.hpp"
#include "Simulator.hpp"

#ifdef SYCL
#include "SyclAllocator.hpp"
#include <sycl/sycl.hpp>
#else
#include "HostAllocator.hpp"
#endif

#include <mpi.h>
#include <Parallel/MPI.h>

#include <fstream>
#include <memory>
#include <iostream>
#include <iterator>
#include <sstream>

int main(int argc, char** argv) {
  auto args = seissol::Args{};
  try {
    args = seissol::ArgParser::parseArgs(argc, argv);
    if (args.showHelp || args.prefix.empty()) {
      seissol::ArgParser::showHelp(std::cout);
      return 0;
    }
  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  seissol::MPI::mpi.init(argc, argv);
  const int rank = seissol::MPI::mpi.rank();

  auto oss = std::ostringstream{};
  oss << args.prefix << "." << rank << ".json";

  auto f = std::ifstream(oss.str());
  if (!f.is_open()) {
    std::cerr << "Could not open " << oss.str() << " for reading." << std::endl;
    MPI_Abort(seissol::MPI::mpi.comm(), -1);
    return -1;
  }

  auto json = std::string(std::istreambuf_iterator<char>{f}, {});
  f.close();
  auto clusters = seissol::LtsStructureParser(json, &std::cerr).parse();
  if (!clusters) {
    std::cerr << "Could not parse " << oss.str() << std::endl;
    MPI_Abort(seissol::MPI::mpi.comm(), -1);
    return -1;
  }

  auto const size = seissol::MPI::mpi.size();
  for (auto const& tc : *clusters) {
    for (auto const& reg : tc.regions) {
      if (reg.neighborRank >= size) {
        std::cerr << "Detected region with neighbor rank " << reg.neighborRank
                  << " but the comm size is only " << size << std::endl;
        std::cerr << "Did you start the comm proxy with the same number of ranks as you used for"
                  << std::endl
                  << "the LTS structure dump?" << std::endl;
        MPI_Abort(seissol::MPI::mpi.comm(), -1);
        return -1;
      }
    }
  }

  if (rank == 0) {
    std::cout << "====== SeisSol communication reproducer ======" << std::endl;
    std::cout << "JSON dump:            " << args.prefix << std::endl;
    std::cout << "Memory type:          " << to_string(args.memType) << std::endl;
    std::cout << "Synchronization time: " << args.synchronizationTime << std::endl;
  }

#ifdef SYCL
  auto alloc = std::make_shared<seissol::SyclAllocator>(args.memType, sycl::queue{});
#else
  auto alloc = std::make_shared<seissol::HostAllocator>();
#endif

  auto simulator = seissol::Simulator(*clusters, std::move(alloc));
  simulator.simulate(args.synchronizationTime);

  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());
  if (rank == 0) {
    std::cout << "Successfully reached sync time (" << args.synchronizationTime << " s)"
              << std::endl;
  }

  seissol::MPI::mpi.finalize();

  return 0;
}
