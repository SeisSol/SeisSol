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

#ifndef MEMORYMANAGER_H_
#define MEMORYMANAGER_H_

#include <Common/configs.hpp>
#include <Common/templating.hpp>
#include <Initializer/Layout/Memory.hpp>
#include <Initializer/tree/LTSForest.hpp>
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <utils/logger.h>

#include <Initializer/typedefs.hpp>
#include "MemoryAllocator.h"

#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Lut.hpp>
#include <Initializer/DynamicRupture.h>
#include <Initializer/InputAux.hpp>
#include <Initializer/Boundary.h>
#include <Initializer/ParameterDB.h>
#include <Initializer/time_stepping/LtsParameters.h>

#include <Physics/InitialField.h>

#include <vector>
#include <memory>

#include <DynamicRupture/Factory.h>
#include <yaml-cpp/yaml.h>

namespace seissol {
  namespace initializers {
    class MemoryManager {
  private: // explicit private for unit tests
    //! memory allocator
    seissol::memory::ManagedAllocator m_memoryAllocator;

    seissol::initializer::MemoryContainer memoryContainer;

    std::vector<std::unique_ptr<physics::InitialField>> m_iniConds;

    template<typename Config>
    using ConfigFrictionSolver = std::unique_ptr<seissol::dr::friction_law::FrictionSolver<Config>>;
    template<typename Config>
    using ConfigFaultOutputManager = std::unique_ptr<seissol::dr::output::OutputManager<Config>>;

    TransformVariadicT<ConfigFrictionSolver, SupportedConfigs> m_FrictionLaw = nullptr;
    TransformVariadicT<ConfigFaultOutputManager, SupportedConfigs> m_faultOutputManager = nullptr;
    std::shared_ptr<dr::DRParameters> m_dynRupParameters = nullptr;
    std::shared_ptr<YAML::Node> m_inputParams = nullptr;
    std::shared_ptr<time_stepping::LtsParameters> ltsParameters = nullptr;

    EasiBoundary m_easiBoundary;


  public:
    /**
     * Constructor
     **/
    MemoryManager() {}

    /**
     * Destructor, memory is freed by managed allocator
     **/
    ~MemoryManager() {}

    inline seissol::initializer::MemoryContainer& getMemoryContainer() {
      return memoryContainer;
    }

    inline const seissol::initializer::MemoryContainer& getMemoryContainer() const {
      return memoryContainer;
    }

    inline void setInitialConditions(std::vector<std::unique_ptr<physics::InitialField>>&& iniConds) {
      m_iniConds = std::move(iniConds);
    }

    inline const std::vector<std::unique_ptr<physics::InitialField>>& getInitialConditions() {
      return m_iniConds;
    }

    void initializeEasiBoundaryReader(const char* fileName);

    inline EasiBoundary* getEasiBoundaryReader() {
      return &m_easiBoundary;
    }

    inline seissol::dr::friction_law::FrictionSolver* getFrictionLaw() {
        return m_FrictionLaw.get();
    }
    inline seissol::dr::output::OutputManager* getFaultOutputManager() {
        return m_faultOutputManager.get();
    }
    inline seissol::dr::DRParameters* getDRParameters() {
        return m_dynRupParameters.get();
    }

    inline time_stepping::LtsParameters* getLtsParameters() {
        return ltsParameters.get();
    };

    void setInputParams(std::shared_ptr<YAML::Node> params) {
      m_inputParams = params;
      m_dynRupParameters = dr::readParametersFromYaml(m_inputParams);
      ltsParameters = std::make_shared<time_stepping::LtsParameters>(time_stepping::readLtsParametersFromYaml(m_inputParams));
    }

    std::string getOutputPrefix() const {
      return getUnsafe<std::string>((*m_inputParams)["output"], "outputfile");
    }

    bool isLoopStatisticsNetcdfOutputOn() const {
      return getWithDefault((*m_inputParams)["output"], "loopstatisticsnetcdfoutput", false);
    }

#ifdef ACL_DEVICE
  void recordExecutionPaths(bool usePlasticity);
#endif

  void initializeFrictionLaw();
  void initFaultOutputManager();
  void initFrictionData();
};
    bool isAcousticSideOfElasticAcousticInterface(CellMaterialData &material,
                                                  unsigned int face);
    bool isElasticSideOfElasticAcousticInterface(CellMaterialData &material,
                                                 unsigned int face);
    bool isAtElasticAcousticInterface(CellMaterialData &material, unsigned int face);

    bool requiresDisplacement(CellLocalInformation cellLocalInformation,
                              CellMaterialData &material,
                              unsigned int face);
    bool requiresNodalFlux(FaceType f);
    }
}

#endif
