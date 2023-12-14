// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SIMULATOR_20231020_HPP
#define SIMULATOR_20231020_HPP

#include "Allocator.hpp"
#include "LtsStructure.hpp"

#include <Initializer/MeshStructure.hpp>
#include <Solver/time_stepping/AbstractTimeCluster.h>
#include <Solver/time_stepping/AbstractGhostTimeCluster.h>

#include <memory>
#include <vector>

namespace seissol {

class Simulator {
  public:
  Simulator(std::vector<TimeCluster> const& ltsStructure, std::shared_ptr<Allocator> alloc);
  ~Simulator();

  Simulator(Simulator const&) = delete;
  Simulator(Simulator&&) = delete;
  Simulator& operator=(Simulator const&) = delete;
  Simulator& operator=(Simulator&&) = delete;

  void simulate(double synchronizationTime);

  private:
  bool poll();

  std::shared_ptr<Allocator> alloc_;
  std::vector<MeshStructure> meshStructure_;
  real* buffer_;
  std::vector<std::unique_ptr<time_stepping::AbstractTimeCluster>> clusters_;
  std::vector<std::unique_ptr<time_stepping::AbstractGhostTimeCluster>> ghostClusters_;
};

} // namespace seissol

#endif // SIMULATOR_20231020_HPP
