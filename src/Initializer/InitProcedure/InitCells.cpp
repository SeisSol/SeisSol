
#include <vector>
#include "Initializer/ParameterDB.h"
#include "Initializer/InputParameters.hpp"
#include "Initializer/CellLocalMatrices.h"
#include "Initializer/LTS.h"
#include "Physics/Attenuation.hpp"
#include "Equations/datastructures.hpp"

#include "SeisSol.h"
#include "Init.hpp"
#include "InitCells.hpp"

#include <cmath>

using namespace seissol::initializer;

using Material = seissol::model::MaterialClass;
using Plasticity = seissol::model::Plasticity;

template<typename T>
std::vector<T> queryDB(seissol::initializers::QueryGenerator* queryGen, const std::string& fileName) {
  const auto& meshReader = seissol::SeisSol::main.meshReader();
    std::vector<T> vectorDB(meshReader.getElements().size());
    seissol::initializers::MaterialParameterDB<T> parameterDB;
    parameterDB.setMaterialVector(&vectorDB);
    parameterDB.evaluateModel(fileName, queryGen);
    return vectorDB;
}

void initializeCellMaterial(seissol::initializer::initprocedure::LtsInfo& ltsInfo) {
    const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();
    const auto& meshReader = seissol::SeisSol::main.meshReader();
    initializers::MemoryManager& memoryManager = seissol::SeisSol::main.getMemoryManager();

    seissol::initializers::QueryGenerator* queryGen = seissol::initializers::getBestQueryGenerator(
      seissol::initializer::parameters::modelAnelastic(),
      ssp.model.plasticity,
      seissol::initializer::parameters::modelAnisotropic(),
      seissol::initializer::parameters::modelPoroelastic(),
      ssp.model.useCellHomogenizedMaterial,
      meshReader);
    auto materialsDB = queryDB<Material>(queryGen, ssp.model.materialFileName);
    std::vector<Plasticity> plasticityDB;
    if (ssp.model.plasticity) {
        plasticityDB = queryDB<Plasticity>(queryGen, ssp.model.materialFileName);
    }

#if defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
    // we need to compute all model parameters before we can use them...
    // TODO: integrate this with the Viscoelastic material class or the ParameterDB directly?
    logInfo() << "Initializing attenuation.";
#ifdef OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < materialsDB.size(); ++i) {
        auto& cellmat = materialsDB[i];
        seissol::physics::fitAttenuation(cellmat, ssp.model.freqCentral, ssp.model.freqRatio);
    }
#endif

    logInfo() << "Setting cell materials in the LTS tree.";
#ifdef OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < meshReader.getElements().size(); ++i) {
        auto& cellInformation = memoryManager.getLtsLut()->lookup(memoryManager.getLts()->cellInformation, i);
        auto& material = memoryManager.getLtsLut()->lookup(memoryManager.getLts()->material, i);

        // set the materials for the cell volume and its faces
        material.local = materialsDB[i];
        for (size_t side = 0; side < 4; ++side) {
            if (cellInformation.faceTypes[side] == FaceType::regular || cellInformation.faceTypes[side] == FaceType::periodic || cellInformation.faceTypes[side] == FaceType::dynamicRupture) {
                // use the neighbor face material info in case that we are not at a boundary
                auto neighbor = cellInformation.faceNeighborIds[side];
                material.neighbor[side] = materialsDB[neighbor]; // TODO: does this work?
            }
            else {
                // otherwise, use the material from the own cell
                material.neighbor[side] = material.local;
            }
        }

        // if enabled, set up the plasticity as well
        if (ssp.model.plasticity) {
            auto& plasticity = memoryManager.getLtsLut()->lookup(memoryManager.getLts()->plasticity, i);

            plasticity.initialLoading[0] = plasticityDB[i].s_xx;
            plasticity.initialLoading[1] = plasticityDB[i].s_yy;
            plasticity.initialLoading[2] = plasticityDB[i].s_zz;
            plasticity.initialLoading[3] = plasticityDB[i].s_xy;
            plasticity.initialLoading[4] = plasticityDB[i].s_yz;
            plasticity.initialLoading[5] = plasticityDB[i].s_xz;

            double angularFriction = std::atan(plasticityDB[i].bulkFriction);

            plasticity.cohesionTimesCosAngularFriction = plasticityDB[i].plastCo * std::cos(angularFriction);
            plasticity.sinAngularFriction = std::sin(angularFriction);
#ifndef USE_ANISOTROPIC
            plasticity.mufactor = 1.0 / (2.0 * material.local.mu);
#else
            plasticity.mufactor = 3.0 / (2.0 * (material.local.c44 + material.local.c55 + material.local.c66));
#endif
        }
    }

    // set tv for all time clusters
    if (ssp.model.plasticity) {
      seissol::SeisSol::main.timeManager().setTv(ssp.model.tv);
    }
}

void initializeCellMatrices(seissol::initializer::initprocedure::LtsInfo& ltsInfo)
{
  const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();

  // \todo Move this to some common initialization place
  auto& meshReader = seissol::SeisSol::main.meshReader();
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();

  seissol::initializers::initializeCellLocalMatrices( meshReader,
                                                      memoryManager.getLtsTree(),
                                                      memoryManager.getLts(),
                                                      memoryManager.getLtsLut(),
                                                      ltsInfo.timeStepping);

  seissol::initializers::initializeDynamicRuptureMatrices( meshReader,
                                                           memoryManager.getLtsTree(),
                                                           memoryManager.getLts(),
                                                           memoryManager.getLtsLut(),
                                                           memoryManager.getDynamicRuptureTree(),
                                                           memoryManager.getDynamicRupture(),
                                                           ltsInfo.ltsMeshToFace,
                                                           *memoryManager.getGlobalDataOnHost(),
                                                           ltsInfo.timeStepping );

  memoryManager.initFrictionData();
  seissol::SeisSol::main.getMemoryManager().getFaultOutputManager()->initFaceToLtsMap();

  seissol::initializers::initializeBoundaryMappings(meshReader,
                                                    memoryManager.getEasiBoundaryReader(),
                                                    memoryManager.getLtsTree(),
                                                    memoryManager.getLts(),
                                                    memoryManager.getLtsLut());

#ifdef ACL_DEVICE
  initializers::copyCellMatricesToDevice(memoryManager.getLtsTree(),
                                         memoryManager.getLts(),
                                         memoryManager.getDynamicRuptureTree(),
                                         memoryManager.getDynamicRupture(),
                                         memoryManager.getBoundaryTree(),
                                         memoryManager.getBoundary());

  memoryManager.recordExecutionPaths(ssp.model.plasticity);
#endif
}

template<typename T>
void synchronize(seissol::initializer::initprocedure::LtsInfo& ltsInfo, seissol::initializers::Variable<T> const& handle)
{
  initializers::MemoryManager& memoryManager = seissol::SeisSol::main.getMemoryManager();
  unsigned *const (&meshToLts)[seissol::initializers::Lut::MaxDuplicates] = memoryManager.getLtsLut()->getMeshToLtsLut(handle.mask);
  unsigned* duplicatedMeshIds = memoryManager.getLtsLut()->getDuplicatedMeshIds(handle.mask);
  unsigned numberOfDuplicatedMeshIds = memoryManager.getLtsLut()->getNumberOfDuplicatedMeshIds(handle.mask);
  T* var = memoryManager.getLtsTree()->var(handle);
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned dupMeshId = 0; dupMeshId < numberOfDuplicatedMeshIds; ++dupMeshId) {
    unsigned meshId = duplicatedMeshIds[dupMeshId];
    T* ref = &var[ meshToLts[0][meshId] ];
    for (unsigned dup = 1; dup < seissol::initializers::Lut::MaxDuplicates && meshToLts[dup][meshId] != std::numeric_limits<unsigned>::max(); ++dup) {
      memcpy(reinterpret_cast<void*>(&var[ meshToLts[dup][meshId] ]),
             reinterpret_cast<void*>(ref),
             sizeof(T));
    }
  }
}

void synchronizeCellLocalData(seissol::initializer::initprocedure::LtsInfo& ltsInfo) {
  const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();
  initializers::MemoryManager& memoryManager = seissol::SeisSol::main.getMemoryManager();
  const auto& lts = memoryManager.getLts();
  synchronize(ltsInfo, lts->material);
  if (ssp.model.plasticity) {
    synchronize(ltsInfo, lts->plasticity);
  }
  synchronize(ltsInfo, lts->dofs);
  if (kernels::size<tensor::Qane>() > 0) {
    synchronize(ltsInfo, lts->dofsAne);
  }
}

void seissol::initializer::initprocedure::initCells(seissol::initializer::initprocedure::LtsInfo& ltsInfo) {
    SCOREP_USER_REGION("init_model", SCOREP_USER_REGION_TYPE_FUNCTION);

    logInfo() << "Begin init model.";

    // Call the pre mesh initialization hook
	seissol::Modules::callHook<seissol::PRE_MODEL>();

    seissol::Stopwatch watch;
	watch.start();
  
  logInfo() << "Initialize cell material parameters.";
	initializeCellMaterial(ltsInfo);
  logInfo() << "Initialize cell-local matrices.";
  initializeCellMatrices(ltsInfo);
  logInfo() << "Synchronize data structures.";
  synchronizeCellLocalData(ltsInfo);

	watch.pause();
	watch.printTime("Model initialized in:");

    // Call the post mesh initialization hook
	seissol::Modules::callHook<seissol::POST_MODEL>();

    logInfo() << "End init model.";
}
