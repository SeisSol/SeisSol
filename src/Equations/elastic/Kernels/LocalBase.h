/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Volume kernel of SeisSol.
 **/

#ifndef KERNELS_VOLUMEBASE_H_
#define KERNELS_VOLUMEBASE_H_

#include <memory>
#include "generated_code/kernel.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop
#include "Physics/InitialField.h"

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol {
  namespace kernels {
    class LocalBase;
  }
}
struct GlobalData;

class seissol::kernels::LocalBase {
  protected:
    static void checkGlobalData(GlobalData const* global, size_t alignment);
    kernel::volume m_volumeKernelPrototype;
    kernel::localFlux m_localFluxKernelPrototype;
    kernel::localFluxNodal m_nodalLfKrnlPrototype;

    kernel::projectToNodalBoundary m_projectKrnlPrototype;
    kernel::projectToNodalBoundaryRotated m_projectRotatedKrnlPrototype;

    kernels::DirichletBoundary dirichletBoundary;

#ifdef ACL_DEVICE
    kernel::gpu_volume deviceVolumeKernelPrototype;
    kernel::gpu_localFlux deviceLocalFluxKernelPrototype;
    kernel::gpu_localFluxNodal deviceNodalLfKrnlPrototype;
    kernel::gpu_projectToNodalBoundaryRotated deviceProjectRotatedKrnlPrototype;
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

    const std::vector<std::unique_ptr<physics::InitialField>> *initConds;
public:
    virtual void setInitConds(decltype(initConds) initConds) {
      this->initConds = initConds;
    }

    physics::InitialField* getInitCond(size_t index) {
      const auto& condition = this->initConds->at(index);
      return condition.get();
    }
};
#endif

