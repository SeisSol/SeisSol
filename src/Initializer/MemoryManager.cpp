/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Memory management of SeisSol.
 **/
#include "SeisSol.h"
#include "MemoryManager.h"
#include "GlobalData.h"
#include <yateto.h>

#include <Kernels/common.hpp>
#include <generated_code/tensor.h>
#include <unordered_set>
#include <cmath>
#include <type_traits>

#include "tree/LTSForest.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ACL_DEVICE
#include "BatchRecorders/Recorders.h"
#include "device.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#endif // ACL_DEVICE

void seissol::initializers::MemoryManager::initializeEasiBoundaryReader(const char* fileName) {
  const auto fileNameStr = std::string{fileName};
  if (fileNameStr != "") {
    m_easiBoundary = EasiBoundary(fileNameStr);
  }
}


#ifdef ACL_DEVICE
void seissol::initializers::MemoryManager::recordExecutionPaths(bool usePlasticity) {
  recording::CompositeRecorder<seissol::initializers::LTS> recorder;
  recorder.addRecorder(new recording::LocalIntegrationRecorder);
  recorder.addRecorder(new recording::NeighIntegrationRecorder);

  if (usePlasticity) {
    recorder.addRecorder(new recording::PlasticityRecorder);
  }

  container.cluster.visit([&](auto&& ltsview) {
    for (ltsview.tree::leaf_iterator it = ltsview.tree.beginLeaf(Ghost); it != ltsview.tree.endLeaf(); ++it) {
      recorder.record(lts, *it);
    }
  });

  recording::CompositeRecorder<seissol::initializers::DynamicRupture> drRecorder;
  drRecorder.addRecorder(new recording::DynamicRuptureRecorder);
  container.dynrup.visit([&](auto&& dynrupview) {
    for (ltsview.tree::leaf_iterator it = dynrupview.tree.beginLeaf(Ghost); it != dynrupview.tree.endLeaf(); ++it) {
      drRecorder.record(*dynrupview.lts, *it);
    }
  });
}
#endif // ACL_DEVICE

bool seissol::initializers::isAcousticSideOfElasticAcousticInterface(CellMaterialData &material,
                                              unsigned int face) {
#ifdef USE_ANISOTROPIC
  return false;
#else
  constexpr auto eps = std::numeric_limits<real>::epsilon();
  return material.neighbor[face]->getMu() > eps && material.local->getMu() < eps;
#endif
}
bool seissol::initializers::isElasticSideOfElasticAcousticInterface(CellMaterialData &material,
                                             unsigned int face) {
#ifdef USE_ANISOTROPIC
  return false;
#else
  constexpr auto eps = std::numeric_limits<real>::epsilon();
  return material.local->getMu() > eps && material.neighbor[face]->getMu() < eps;
#endif
}

bool seissol::initializers::isAtElasticAcousticInterface(CellMaterialData &material, unsigned int face) {
  // We define the interface cells as all cells that are in the elastic domain but have a
  // neighbor with acoustic material.
#ifndef USE_ANISOTROPIC
  return isAcousticSideOfElasticAcousticInterface(material, face) || isElasticSideOfElasticAcousticInterface(material, face);
#else
  return false;
#endif
}


bool seissol::initializers::requiresDisplacement(CellLocalInformation cellLocalInformation,
                                                 CellMaterialData &material,
                                                 unsigned int face) {
  const auto faceType = cellLocalInformation.faceTypes[face];
  return faceType == FaceType::freeSurface
  || faceType == FaceType::freeSurfaceGravity
  || isAtElasticAcousticInterface(material, face);
}

bool seissol::initializers::requiresNodalFlux(FaceType f) {
  return (f == FaceType::freeSurfaceGravity
          || f == FaceType::dirichlet
          || f == FaceType::analytical);
}

void seissol::initializers::MemoryManager::initializeFrictionLaw() {
  const int rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Initialize Friction Model";

  const auto factory = seissol::dr::factory::getFactory(m_dynRupParameters);
  auto product = factory->produce();
  m_dynRup = std::move(product.ltsview.tree);
  m_DRInitializer = std::move(product.initializer);
  m_FrictionLaw = std::move(product.frictionLaw);
  m_faultOutputManager = std::move(product.output);
}

void seissol::initializers::MemoryManager::initFaultOutputManager() {
  // TODO: switch m_dynRup to shared or weak pointer
  if (m_dynRupParameters->isDynamicRuptureEnabled) {
    /*m_faultOutputManager->setInputParam(*m_inputParams, seissol::SeisSol::main.meshReader());
    m_faultOutputManager->setLtsData(&m_clusterforest,
                                     &m_ltsLut,
                                     &m_dynrupforest);
    m_faultOutputManager->init();*/

  }
}


void seissol::initializers::MemoryManager::initFrictionData() {
  if (m_dynRupParameters->isDynamicRuptureEnabled) {

    m_DRInitializer->initializeFault(m_dynrupforest);

#ifdef ACL_DEVICE
    if (auto* impl = dynamic_cast<dr::friction_law::gpu::FrictionSolverInterface*>(m_FrictionLaw.get())) {
      impl->initSyclQueue();

      LayerMask mask = seissol::initializers::LayerMask(Ghost);
      auto maxSize = m_dynrupview.tree.getMaxClusterSize(mask);
      impl->setMaxClusterSize(maxSize);

      impl->allocateAuxiliaryMemory();
      impl->copyStaticDataToDevice();
    }
#endif // ACL_DEVICE
  }
}

