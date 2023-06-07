#pragma once

#include "WavePropagation/dispatcher.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Kernels/Time.h"
#include "Kernels/Local.h"
#include "Kernels/Neighbor.h"
#include "SeisSol.h"

#include <vector>

namespace seissol::waveprop {
    // template<bool Plasticity> <-- one day, templates may be added here
    class WavePropDispatcherPre : public WavePropDispatcherBase {
    public:
        WavePropDispatcherPre(const seissol::initializers::LTS& lts, seissol::initializers::Layer& layer)
            : lts(lts), layer(layer), globalData(seissol::SeisSol::main.getMemoryManager().getGlobalData())
        {
          timeKernel.setGlobalData(globalData);
          localKernel.setGlobalData(globalData);
          localKernel.setInitConds(&seissol::SeisSol::main.getMemoryManager().getInitialConditions());
          neighborKernel.setGlobalData(globalData);
        }

        virtual void computePredictFlops(long long int& flopsNonZero, long long int& flopsHardware) override {
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
                    GravitationalFreeSurfaceBc::getFlopsDisplacementFace(face,
                                                                        cellInformation[cell].faceTypes[face]);
                    flopsNonZero += nonZeroFlopsDisplacement;
                    flopsHardware += hardwareFlopsDisplacement;
                }
                }
            }
        }

        virtual void computeCorrectFlops(long long int& flopsNonZero, long long int& flopsHardware, long long int& drFlopsNonZero, long long int& drFlopsHardware) override {
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
                                                        cellDRHardware );
                flopsNonZero += cellNonZero;
                flopsHardware += cellHardware;
                drFlopsNonZero += cellDRNonZero;
                drFlopsHardware += cellDRHardware;

                /// \todo add lts time integration
                /// \todo add plasticity
            }
        }

        virtual void setTV(double tv) override {
            this->tv = tv;
        }

    protected: // for now, at least
        CompoundGlobalData globalData;
        seissol::initializers::Layer& layer;
        const seissol::initializers::LTS& lts;
        kernels::Time timeKernel;
        kernels::Local localKernel;
        kernels::Neighbor neighborKernel;
        double tv;
    };
}
