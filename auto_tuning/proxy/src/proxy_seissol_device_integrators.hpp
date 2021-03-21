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

#include <generated_code/tensor.h>
#include <device.h>
#include <stdio.h>

namespace tensor = seissol::tensor;
namespace kernels = seissol::kernels;

namespace proxy::device {
  using deviceT = ::device::DeviceInstance;
  void computeAderIntegration() {
    const deviceT &device = deviceT::getInstance();
    auto& layer = m_ltsTree->child(0).child<Interior>();

    kernels::LocalData::Loader loader;
    loader.load(m_lts, layer);
    kernels::LocalTmp tmp;

    ConditionalBatchTableT &table = layer.getCondBatchTable();

    m_timeKernel.computeBatchedAder(static_cast<double>(m_timeStepWidthSimulation), tmp, table);
    device.api->synchDevice();
  }

  void computeLocalWithoutAderIntegration() {
    const deviceT &device = deviceT::getInstance();
    auto& layer = m_ltsTree->child(0).child<Interior>();
    kernels::LocalData::Loader loader;
    loader.load(m_lts, layer);
    kernels::LocalTmp tmp;

    ConditionalBatchTableT &table = layer.getCondBatchTable();

    m_localKernel.computeBatchedIntegral(table, tmp);
    device.api->synchDevice();
  }

  void computeLocalIntegration() {
    const deviceT &device = deviceT::getInstance();
    auto& layer = m_ltsTree->child(0).child<Interior>();

    kernels::LocalData::Loader loader;
    loader.load(m_lts, layer);
    kernels::LocalTmp tmp;

    ConditionalBatchTableT &table = layer.getCondBatchTable();

    m_timeKernel.computeBatchedAder(static_cast<double>(m_timeStepWidthSimulation), tmp, table);
    m_localKernel.computeBatchedIntegral(table, tmp);
    device.api->synchDevice();
  }

  void computeNeighboringIntegration() {
    const deviceT &device = deviceT::getInstance();
    auto& layer = m_ltsTree->child(0).child<Interior>();

    kernels::NeighborData::Loader loader;
    loader.load(m_lts, layer);

    ConditionalBatchTableT &table = layer.getCondBatchTable();

    seissol::kernels::TimeCommon::computeBatchedIntegrals(m_timeKernel,
                                                          0.0,
                                                         static_cast<double>(m_timeStepWidthSimulation),
                                                         table);
    m_neighborKernel.computeBatchedNeighborsIntegral(table);
    device.api->synchDevice();
  }

  void computeDynRupGodunovState() {
    auto& layer = m_dynRupTree->child(0).child<Interior>();

    ConditionalBatchTableT &table = layer.getCondBatchTable();
    m_dynRupKernel.batchedSpaceTimeInterpolation(table);
  }
} // namespace proxy::device


#endif //SEISSOL_PROXY_SEISSOL_DEVICE_INTEGRATORS_HPP
