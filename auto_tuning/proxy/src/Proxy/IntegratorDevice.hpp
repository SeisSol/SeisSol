/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef SEISSOL_PROXY_SEISSOL_DEVICE_INTEGRATORS_HPP
#define SEISSOL_PROXY_SEISSOL_DEVICE_INTEGRATORS_HPP

#include "generated_code/tensor.h"
#include <device.h>
#include <stdio.h>

namespace tensor = seissol::tensor;
namespace kernels = seissol::kernels;

namespace proxy::device {
using deviceType = ::device::DeviceInstance;
void computeAderIntegration() {
  auto& layer = ltsTree->child(0).child<Interior>();

  kernels::LocalData::Loader loader;
  loader.load(lts, layer);
  kernels::LocalTmp tmp(9.81);

  auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.getConditionalTable<inner_keys::Material>();

  const double timeStepWidth = static_cast<double>(seissol::miniSeisSolTimeStep);
  ComputeGraphType graphType{ComputeGraphType::LocalIntegral};
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, false);

  runtime->runGraph(computeGraphKey, layer, [&](auto& runtime) {
    timeKernel.computeBatchedAder(timeStepWidth, tmp, dataTable, materialTable, false, runtime);
  });
}

void computeLocalWithoutAderIntegration() {
  auto& layer = ltsTree->child(0).child<Interior>();
  kernels::LocalData::Loader loader;
  loader.load(lts, layer);
  kernels::LocalTmp tmp(9.81);

  auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.getConditionalTable<inner_keys::Material>();
  auto& indicesTable = layer.getConditionalTable<inner_keys::Indices>();

  const double timeStepWidth = 0.0;
  ComputeGraphType graphType{ComputeGraphType::LocalIntegral};
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, false);

  runtime->runGraph(computeGraphKey, layer, [&](auto& runtime) {
    localKernel.computeBatchedIntegral(
        dataTable, materialTable, indicesTable, loader, tmp, timeStepWidth, runtime);
  });
}

void computeLocalIntegration() {
  auto& layer = ltsTree->child(0).child<Interior>();

  kernels::LocalData::Loader loader;
  loader.load(lts, layer);
  kernels::LocalTmp tmp(9.81);

  auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.getConditionalTable<inner_keys::Material>();
  auto& indicesTable = layer.getConditionalTable<inner_keys::Indices>();

  const double timeStepWidth = static_cast<double>(seissol::miniSeisSolTimeStep);
  ComputeGraphType graphType{ComputeGraphType::LocalIntegral};
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, false);
  runtime->runGraph(computeGraphKey, layer, [&](auto& runtime) {
    timeKernel.computeBatchedAder(timeStepWidth, tmp, dataTable, materialTable, false, runtime);
    localKernel.computeBatchedIntegral(
        dataTable, materialTable, indicesTable, loader, tmp, 0.0, runtime);
  });
}

void computeNeighboringIntegration() {
  auto& layer = ltsTree->child(0).child<Interior>();

  kernels::NeighborData::Loader loader;
  loader.load(lts, layer);

  const double timeStepWidth = static_cast<double>(seissol::miniSeisSolTimeStep);
  auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();

  seissol::kernels::TimeCommon::computeBatchedIntegrals(
      timeKernel, 0.0, timeStepWidth, dataTable, *runtime);

  ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
  auto computeGraphKey = initializer::GraphKey(graphType);
  runtime->runGraph(computeGraphKey, layer, [&](auto& runtime) {
    neighborKernel.computeBatchedNeighborsIntegral(dataTable, runtime);
  });
}

void computeDynRupGodunovState() {
  auto& layer = dynRupTree->child(0).child<Interior>();

  auto& dataTable = layer.getConditionalTable<inner_keys::Dr>();

  ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
  auto computeGraphKey = initializer::GraphKey(graphType, 0.0);
  runtime->runGraph(computeGraphKey, layer, [&](auto& runtime) {
    dynRupKernel.batchedSpaceTimeInterpolation(dataTable, runtime);
  });
}
} // namespace proxy::device

#endif // SEISSOL_PROXY_SEISSOL_DEVICE_INTEGRATORS_HPP
