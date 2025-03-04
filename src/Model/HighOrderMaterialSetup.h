#pragma once

#include "Model/Common.h"
#include "Model/CommonDatastructures.h"
#include "Model/HighOrderMaterial.h"
#include "generated_code/init.h"
#include <Equations/Datastructures.h>
#include <Equations/acoustic/Model/Datastructures.h>
#include <Equations/acoustic/Model/IntegrationData.h>
#include <generated_code/kernel.h>
#include <generated_code/tensor.h>

#include <utils/logger.h>

namespace seissol::model {
template <typename BaseMaterialT, std::size_t Order>
struct MaterialSetup<HighOrderMaterial<BaseMaterialT, Order>> {
  using MaterialT = HighOrderMaterial<BaseMaterialT, Order>;

  template <typename T>
  static void getTransposedCoefficientMatrix(const MaterialT& material, unsigned dim, T& matM) {
    real starMData[seissol::tensor::starM::size()]{};
    auto starM = seissol::init::starM::view::create(starMData);
    for (std::size_t i = 0; i < MaterialT::Samples3D; ++i) {
      auto subM = starM.subtensor(i, yateto::slice<>(), yateto::slice<>());
      ::seissol::model::getTransposedCoefficientMatrix<BaseMaterialT>(
          material.materials[i], dim, subM);
    }
    seissol::kernel::homStar project;
    project.homproject = seissol::init::homproject::Values;
    project.starM = starMData;
    project.star(0) = matM.data();
    project.execute();
  }

  template <typename Tloc, typename Tneigh>
  static void getTransposedGodunovState(const MaterialT& local,
                                        const MaterialT& neighbor,
                                        FaceType faceType,
                                        Tloc& qGodLocal,
                                        Tneigh& qGodNeighbor) {
    real godLocalMData[seissol::tensor::godLocalM::size()]{};
    real godNeighborMData[seissol::tensor::godNeighborM::size()]{};
    auto godLocalM = seissol::init::godLocalM::view::create(godLocalMData);
    auto godNeighborM = seissol::init::godNeighborM::view::create(godNeighborMData);
    for (std::size_t i = 0; i < MaterialT::Samples3D; ++i) {
      auto subLocal = godLocalM.subtensor(i, yateto::slice<>(), yateto::slice<>());
      auto subNeighbor = godNeighborM.subtensor(i, yateto::slice<>(), yateto::slice<>());
      ::seissol::model::getTransposedGodunovState<BaseMaterialT>(
          local.materials[i], neighbor.materials[i], faceType, subLocal, subNeighbor);
    }
    seissol::kernel::homFluxsolver project;
    project.homproject = seissol::init::homproject::Values;
    project.godLocalM = godLocalMData;
    project.godNeighborM = godNeighborMData;
    project.QgodLocal = qGodLocal.data();
    project.QgodNeighbor = qGodNeighbor.data();
    project.execute();
  }

  static MaterialT getRotatedMaterialCoefficients(real rotationParameters[36],
                                                  MaterialT& material) {
    MaterialT rotated;
    for (std::size_t i = 0; i < MaterialT::Samples3D; ++i) {
      rotated.materials[i] = ::seissol::model::getRotatedMaterialCoefficients(
          rotationParameters, material.materials[i]);
    }
    return rotated;
  }

  // maybe done
  static void initializeSpecificLocalData(const MaterialT& material,
                                          real timeStepWidth,
                                          typename MaterialT::LocalSpecificData* localData) {
    ::seissol::model::initializeSpecificLocalData<BaseMaterialT>(
        material.materials[0], timeStepWidth, localData);
  }

  // maybe done
  static void
      initializeSpecificNeighborData(const MaterialT& material,
                                     typename MaterialT::NeighborSpecificData* neighborData) {
    ::seissol::model::initializeSpecificNeighborData<BaseMaterialT>(material.materials[0],
                                                                    neighborData);
  }

  // maybe done
  static void getPlaneWaveOperator(
      const MaterialT& material,
      const double n[3],
      std::complex<double> mdata[MaterialT::NumQuantities * MaterialT::NumQuantities]) {
    ::seissol::model::getPlaneWaveOperator<BaseMaterialT>(material.materials[0], n, mdata);
  }

  // maybe done
  template <typename T>
  static void getTransposedSourceCoefficientTensor(const MaterialT& material, T& sourceMatrix) {
    ::seissol::model::getTransposedSourceCoefficientTensor<BaseMaterialT>(material.materials[0],
                                                                          sourceMatrix);
  }

  static void getFaceRotationMatrix(const VrtxCoords normal,
                                    const VrtxCoords tangent1,
                                    const VrtxCoords tangent2,
                                    init::T::view::type& matT,
                                    init::Tinv::view::type& matTinv) {
    ::seissol::model::getFaceRotationMatrix<BaseMaterialT>(
        normal, tangent1, tangent2, matT, matTinv);
  }
};
} // namespace seissol::model
