#pragma once

#include "WavePropagation/dispatcher.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Kernels/Time.h"
#include "Kernels/Local.h"
#include "Kernels/Neighbor.h"
#include "SeisSol.h"

#include "Common/cellconfig.hpp"

#include <vector>

namespace seissol::waveprop {

template <typename Config>
class WavePropDispatcherPre : public WavePropDispatcherBase {
  public:
  using RealT = typename Config::RealT;
  using MaterialT = typename Config::MaterialT;

  WavePropDispatcherPre(const seissol::initializers::LTS<Config>& lts,
                        seissol::initializers::Layer& layer)
      : lts(lts), layer(layer),
        globalData(seissol::SeisSol::main.getMemoryManager().getGlobalData()) {
    timeKernel.setGlobalData(globalData);
    localKernel.setGlobalData(globalData);
    localKernel.setInitConds(&seissol::SeisSol::main.getMemoryManager().getInitialConditions());
    neighborKernel.setGlobalData(globalData);
  }

  void computePredictFlops(long long int& flopsNonZero, long long int& flopsHardware) override {
    flopsNonZero = 0;
    flopsHardware = 0;

    auto* cellInformation = layer.var(lts.cellInformation);
    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      unsigned cellNonZero, cellHardware;
      timeKernel.flopsAder(cellNonZero, cellHardware);
      flopsNonZero += cellNonZero;
      flopsHardware += cellHardware;
      localKernel.flopsIntegral(cellInformation[cell].faceTypes, cellNonZero, cellHardware);
      flopsNonZero += cellNonZero;
      flopsHardware += cellHardware;
      // Contribution from displacement/integrated displacement
      for (unsigned face = 0; face < 4; ++face) {
        if (cellInformation->faceTypes[face] == FaceType::freeSurfaceGravity) {
          const auto [nonZeroFlopsDisplacement, hardwareFlopsDisplacement] =
              seissol::kernels::GravitationalFreeSurfaceBc<Config>::getFlopsDisplacementFace(
                  face, cellInformation[cell].faceTypes[face]);
          flopsNonZero += nonZeroFlopsDisplacement;
          flopsHardware += hardwareFlopsDisplacement;
        }
      }
    }
  }

  void computeCorrectFlops(long long int& flopsNonZero,
                           long long int& flopsHardware,
                           long long int& drFlopsNonZero,
                           long long int& drFlopsHardware) override {
    flopsNonZero = 0;
    flopsHardware = 0;
    drFlopsNonZero = 0;
    drFlopsHardware = 0;

    auto* cellInformation = layer.var(lts.cellInformation);
    auto* drMapping = layer.var(lts.drMapping);
    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      unsigned cellNonZero, cellHardware;
      long long cellDRNonZero, cellDRHardware;
      neighborKernel.flopsNeighborsIntegral(cellInformation[cell].faceTypes,
                                            cellInformation[cell].faceRelations,
                                            drMapping[cell],
                                            cellNonZero,
                                            cellHardware,
                                            cellDRNonZero,
                                            cellDRHardware);
      flopsNonZero += cellNonZero;
      flopsHardware += cellHardware;
      drFlopsNonZero += cellDRNonZero;
      drFlopsHardware += cellDRHardware;

      /// \todo add lts time integration
      /// \todo add plasticity
    }
  }

  void setTV(double tv) override { this->tv = tv; }

  protected: // for now, at least
  CompoundGlobalData<Config> globalData;
  seissol::initializers::Layer& layer;
  const seissol::initializers::LTS<Config>& lts;
  seissol::waveprop::kernel::time::Time<Config> timeKernel;
  seissol::waveprop::kernel::local::Local<Config> localKernel;
  seissol::waveprop::kernel::neighbor::Neighbor<Config> neighborKernel;
  double tv;
};
} // namespace seissol::waveprop
