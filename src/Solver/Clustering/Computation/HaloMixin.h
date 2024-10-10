// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <type_traits>
namespace seissol::solver::clustering::computation {

template<typename Base>
class InteriorMixin : public Base {
};

template<typename Base, ComputeStep SendStep>
class CopyMixin : public Base {
  static_assert(std::is_base_of_v<AbstractTimeCluster, Base>, "Base needs to inherit from AbstractTimeCluster");

  void start() override { dataSent = false; }

  void runCompute(ComputeStep step) override {
    if (step == SendStep && dataSent) {
      neighbor->stopTo(this->streamRuntime);
    }

    Base::runCompute(step);

    if (step == SendStep) {
      neighbor->startFrom(this->streamRuntime);
      dataSent = true;
    }
  }

  void reset() override {
    if (dataSent) {
      neighbor->stopTo(this->streamRuntime);
      dataSent = false;
    }
  }
  
  private:
  std::shared_ptr<communication::SendNeighborCluster> neighbor;
  bool dataSent{false};
};

template<ComputeStep SendStep, ComputeStep ArriveStep>
class GhostMixin {

};

}
