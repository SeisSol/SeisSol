#include "Device.hpp"

namespace {

#ifdef ACL_DEVICE
static void seissol::initializers::MemoryManager::deriveRequiredScratchpadMemoryForWp(seissol::initializers::LTSLayerRef<seissol::initializers::LTS<Config>>& layerView) {
  constexpr size_t totalDerivativesSize = yateto::computeFamilySize<tensor::dQ>();
  constexpr size_t nodalDisplacementsSize = tensor::averageNormalDisplacement::size();

  for (auto layer = ltsview.tree.beginLeaf(Ghost); layer != ltsview.tree.endLeaf(); ++layer) {

    CellLocalInformation *cellInformation = layer->var(ltsview.lts.cellInformation);
    std::unordered_set<real *> registry{};
    real *(*faceNeighbors)[4] = layer->var(ltsview.lts.faceNeighbors);

    unsigned derivativesCounter{0};
    unsigned integratedDofsCounter{0};
    unsigned nodalDisplacementsCounter{0};

    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      bool needsScratchMemForDerivatives = (cellInformation[cell].ltsSetup >> 9) % 2 == 0;
      if (needsScratchMemForDerivatives) {
        ++derivativesCounter;
      }
      ++integratedDofsCounter;

      // include data provided by ghost layers
      for (unsigned face = 0; face < 4; ++face) {
        real *neighbourBuffer = faceNeighbors[cell][face];

        // check whether a neighbour element idofs has not been counted twice
        if ((registry.find(neighbourBuffer) == registry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighbourBuffer != nullptr) {
            if (cellInformation[cell].faceTypes[face] != FaceType::outflow &&
                cellInformation[cell].faceTypes[face] != FaceType::dynamicRupture) {

              bool isNeighbProvidesDerivatives = ((cellInformation[cell].ltsSetup >> face) % 2) == 1;
              if (isNeighbProvidesDerivatives) {
                ++integratedDofsCounter;
              }
              registry.insert(neighbourBuffer);
            }
          }
        }

        if (cellInformation[cell].faceTypes[face] == FaceType::freeSurfaceGravity) {
          ++nodalDisplacementsCounter;
        }

      }
    }
    layer->setScratchpadSize(ltsview.lts.integratedDofsScratch,
                             integratedDofsCounter * tensor::I::size() * sizeof(real));
    layer->setScratchpadSize(ltsview.lts.derivativesScratch,
                             derivativesCounter * totalDerivativesSize * sizeof(real));
    layer->setScratchpadSize(ltsview.lts.nodalAvgDisplacements,
                             nodalDisplacementsCounter * nodalDisplacementsSize * sizeof(real));
  }
}

#endif

} // namespace

namespace seissol::initializer::internal {

void allocateScratchpads(MemoryContainer& container) {
#ifdef ACL_DEVICE
    container.cluster.visitLayers([](auto&& layerview) {

    });

    container.cluster.visit([](auto&& treeview) {
        treeview.tree.allocateScratchpadArrays();
    });

    container.cluster.visitLayers([](auto&& layerview) {
        using Config = typename std::decay_t<decltype(layerview)>::ConfigT;
        constexpr size_t IDofsSize = Yateto<Config>::Tensor::Q::size() * sizeof(typename Config::RealT);
        const auto layerSize = layerview.layer.getNumberOfCells();
        layerview.layer.setScratchpadSize(layerview.lts.idofsPlusOnDevice, IDofsSize * layerSize);
        layerview.layer.setScratchpadSize(layerview.lts.idofsMinusOnDevice, IDofsSize * layerSize);
    });

    container.cluster.visit([](auto&& treeview) {
        treeview.tree.allocateScratchpadArrays();
    });
#endif
}

} // namespace seissol::initializer::internal
