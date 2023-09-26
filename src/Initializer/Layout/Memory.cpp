#include "Memory.hpp"


#include "Initializer/Layout/Internal/MeshLayout.hpp"
#include "SeisSol.h"
#include <Common/configs.hpp>
#include <Initializer/CellLocalMatrices.h>
#include <Initializer/ConfigFile.hpp>
#include <Initializer/Layout/Internal/Buckets.hpp>
#include <Initializer/ParameterDB.h>
#include <Initializer/time_stepping/Timestep.hpp>
#include <Initializer/tree/Layer.hpp>
#include <Initializer/tree/LayerMap.hpp>
#include <Initializer/tree/TimeCluster.hpp>
#include <Initializer/typedefs.hpp>
#include <Solver/time_stepping/TimeCluster.h>
#include "Initializer/ParameterMaterialDB.hpp"
#include <cstddef>
#include <functional>
#include <limits>
#include <variant>
#include "MeshLayout.hpp"
#include "Internal/LtsSetup.hpp"

namespace seissol::initializer {

static std::pair<std::vector<int>, std::vector<int>> computeClusters(const geometry::MeshReader& meshReader) {
    double cfl = 0.5; // TODO: import
    double maximumAllowedTimeStep = 5000; // TODO: import
    double rate = 2; // TODO: import
    double wiggle = 1; // TODO: import
    auto globalTimestep = seissol::initializer::computeTimesteps(cfl, maximumAllowedTimeStep, initializers::CellToVertexArray::fromMeshReader(meshReader));
    auto clusters = clusterTimesteps(globalTimestep, rate, wiggle);
    enforceMaximumDifference(1, clusters, meshReader);

    auto ghostClusters = meshReader.exchangeGhostData<int, int>(clusters, MPI::mpi.comm(), MPI_INT);

    return {clusters, ghostClusters};
}

static void initializeCellMatrices(const geometry::MeshReader& meshReader, MemoryContainer& memoryContainer) {
  seissol::initializers::initializeCellLocalMatrices(meshReader,
                                                     memoryContainer.cluster,
                                                     ltsInfo.timeStepping);

  seissol::initializers::initializeDynamicRuptureMatrices(meshReader,
                                                          memoryContainer.cluster,
                                                          memoryContainer.clusterBackmap,
                                                          memoryContainer.dynrup,
                                                          memoryContainer.dynrupBackmap,
                                                          memoryContainer.globalDataStorage,
                                                          ltsInfo.timeStepping);

  memoryManager.initFrictionData();

  seissol::initializers::initializeBoundaryMappings(meshReader,
                                                    memoryManager.getEasiBoundaryReader(),
                                                    memoryContainer.cluster);

#ifdef ACL_DEVICE
  initializers::copyCellMatricesToDevice(memoryManager.cluster,
                                         memoryManager.dynrup,
                                         memoryManager.boundary);

  memoryManager.recordExecutionPaths(seissolParams.model.plasticity);
#endif
}

static void allocateDynRupTree(MemoryContainer& container) {
    container.cluster.visitTwoLayers([&](auto&& ltsview, auto&& dynrupview) {

    }, container.dynrup);

    container.dynrup.allocateTouchVariables();

    container.cluster.visitTwoLayers([&](auto&& ltsview, auto&& dynrupview) {

    }, container.dynrup);
}

static void allocateBoundaryTree(MemoryContainer& container) {
    container.cluster.visitTwoLayers([&](auto&& ltsview, auto&& boundaryview) {
        auto* cellInformation = ltsview.layer.var(ltsview.lts.cellInformation);
        unsigned numberOfBoundaryFaces = 0;
        if (ltsview.icg != Ghost) {
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+ : numberOfBoundaryFaces)
#endif // _OPENMP
            for (unsigned cell = 0; cell < ltsview.layer.getNumberOfCells(); ++cell) {
                for (unsigned face = 0; face < 4; ++face) {
                    if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
                        ++numberOfBoundaryFaces;
                    }
                }
            }
        }
        boundaryview.layer.setNumberOfCells(numberOfBoundaryFaces);
    }, container.boundary);

    container.boundary.allocateTouchVariables();

    container.cluster.visitTwoLayers([&](auto&& ltsview, auto&& boundaryview) {
        auto* cellInformation = ltsview.layer.var(ltsview.lts.cellInformation);
        auto* boundaryMapping = ltsview.layer.var(ltsview.lts.boundaryMapping);
        auto* faceInformation = boundaryview.layer.var(boundaryview.lts.faceInformation);

        std::size_t boundaryFace = 0;
        for (unsigned cell = 0; cell < ltsview.layer.getNumberOfCells(); ++cell) {
            for (unsigned face = 0; face < 4; ++face) {
                if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
                    boundaryMapping[cell][face].nodes = faceInformation[boundaryFace].nodes;
                    boundaryMapping[cell][face].TData = faceInformation[boundaryFace].TData;
                    boundaryMapping[cell][face].TinvData = faceInformation[boundaryFace].TinvData;
                    boundaryMapping[cell][face].easiBoundaryMap = faceInformation[boundaryFace].easiBoundaryMap;
                    boundaryMapping[cell][face].easiBoundaryConstant = faceInformation[boundaryFace].easiBoundaryConstant;
                    ++boundaryFace;
                } else {
                    boundaryMapping[cell][face].nodes = nullptr;
                    boundaryMapping[cell][face].TData = nullptr;
                    boundaryMapping[cell][face].TinvData = nullptr;
                    boundaryMapping[cell][face].easiBoundaryMap = nullptr;
                    boundaryMapping[cell][face].easiBoundaryConstant = nullptr;
                }
            }
        }
    }, container.boundary);
}

MemoryContainer setupMemory(const geometry::MeshReader& meshReader, seissol::time_stepping::TimeManager& timeManager) {
    int rank = MPI::mpi.rank();
    logInfo(rank) << "Starting memory setup...";

    MemoryContainer container;

    auto seissolParameters = seissol::SeisSol::main.getSeisSolParameters();
    auto configs = readConfigFile(seissolParameters.model.configFileName);

    logInfo(rank) << "Compute time clusters...";
    auto [clusters, ghostClusters] = computeClusters(meshReader);
    auto colors = std::vector<int>(meshReader.getElements().size());
    auto ghostColors = std::vector<int>(meshReader.getBoundaryElements().size());

    // TODO: container.colorMap
    const auto clusterCount = seissol::SeisSol::main.maxNumberOfClusters;

    container.colorMap = initializers::ClusterColorMap(initializers::EnumLayer<SupportedConfigs>(configs), initializers::RangeLayer(0, clusterCount - 1));

    for (std::size_t i = 0; i < colors.size(); ++i) {
        colors[i] = container.colorMap.color(configs[meshReader.getElements()[i].group].id, clusters[i]);
    }
    for (std::size_t i = 0; i < ghostColors.size(); ++i) {
        ghostColors[i] = container.colorMap.color(configs[meshReader.getBoundaryElements()[i].get().neighborGroup].id, ghostClusters[i]);
    }

    const auto maxColor = clusterCount * std::variant_size_v<SupportedConfigs>;

    logInfo(rank) << "Compute cell layout...";
    auto meshLayout = seissol::initializer::internal::layoutCells(colors, ghostColors, maxColor, meshReader);

    container.cluster.initialize(clusterCount);
    container.boundary.initialize(clusterCount);
    container.dynrup.initialize(clusterCount);

    // container.clusterBackmap = TODO;
    // container.ghostClusterBackmap = TODO;

    logInfo(rank) << "Setting up cell storage...";
    // pass 1: generate backmap/lut, setup layers, allocate trees
    container.cluster.visitLayers([&](auto&& layerview) {
        const auto& layout = meshLayout[toColor(layerview.config, layerview.cluster)];
        layerview.layer.setNumberOfCells(layout.cells(layerview.icg).size());
    });
    container.cluster.allocateTouchVariables();
    container.cluster.visitLayers([&](auto&& layerview) {
        const auto& layout = meshLayout[toColor(layerview.config, layerview.cluster)];
        const auto& testVariable = layerview.lts.cellInformation; // just need pick a test variable that is available everywhere
        const auto& testTree = layerview.tree.var(testVariable);
        const auto& testLayer = layerview.layer.var(testVariable);
        auto* secondaryCellInformation = layerview.var(layerview.lts.secondaryCellInformation);
        const auto elemSize = testVariable.Size;

        auto addToBackmapCopyInterior = [&](auto cell, auto index){
            return container.clusterBackmap.addElement(layerview.config, testTree, testLayer, cell, index);
        };
        auto addToBackmapGhost = [&](auto cell, auto index){
            return container.ghostClusterBackmap.addElement(layerview.config, testTree, testLayer, cell, index);
        };
        const auto& addToBackmap = ([&](){
            if (layerview.icg == Interior || layerview.icg == Copy) {
                return addToBackmapCopyInterior;
            }
            else {
                return addToBackmapGhost;
            }
        })();

        std::size_t index = 0;
        for (const auto& cell : layout.cells(layerview.icg)) {
            // TODO: two boundary elements with the same source (relevant for Ghost)
            auto dup = addToBackmap(cell, index);

            secondaryCellInformation[index].meshId = cell;
            secondaryCellInformation[index].configId = layerview.config;
            secondaryCellInformation[index].clusterId = layerview.cluster;
            secondaryCellInformation[index].duplicate = dup;
            secondaryCellInformation[index].rank = -1;
            for (int face = 0; face < 4; ++face) {
                secondaryCellInformation[index].neighborRanks[face] = -1;
            }
            ++index;
        }
    });

    // pass 2: setup cellinfo, materials, and LTS setup (first pass)
    

    // material retrieval for interior/copy and ghost layers

    logInfo(rank) << "Retrieving material data for the last time...";
    auto [materialsDB, plasticityDB] = seissol::model::queryMaterial(
        configs, seissol::initializers::CellToVertexArray::fromMeshReader(meshReader), true);
    auto [materialsDBGhost, plasticityDBGhost] = seissol::model::queryMaterial(
        configs,
        seissol::initializers::CellToVertexArray::fromVectors(ghostVertices, ghostGroups),
        true);
    logInfo(rank) << "Setting up cell data...";
    container.cluster.visitLayers([&](auto&& layerview) {
        using Config = typename std::decay_t<decltype(layerview)>::ConfigT;
        using MaterialT = typename Config::MaterialT;

        const auto& layout = meshLayout[toColor(layerview.config, layerview.cluster)];
        auto* materialData = layerview.var(layerview.lts.materialData);
        auto* secondaryCellInformation = layerview.var(layerview.lts.secondaryCellInformation);
        auto* cellInformation = layerview.var(layerview.lts.cellInformation);
        if (layerview.icg == Interior || layerview.icg == Copy) {
            auto* material = layerview.var(layerview.lts.material);
            auto* plasticity = layerview.var(layerview.lts.plasticity);
            std::size_t index = 0;

#ifdef _OPENMP
// #pragma omp parallel for schedule(static)
#endif
            for (const auto& cell : layout.cells(layerview.icg)) {
                new (materialData[index]) MaterialT(materialsDB[cell]); // TODO: initAssign
                new (plasticity[index]) model::PlasticityData<typename Config::RealT>(plasticityDB[cell], &materialData[index]); // TODO: initAssign

                material[index].local = &materialData[index];
                secondaryCellInformation[index].meshId = cell;
                secondaryCellInformation[index].configId = layerview.config;
                secondaryCellInformation[index].clusterId = layerview.cluster;
                secondaryCellInformation[index].rank = MPI::mpi.rank();
                for (int face = 0; face < 4; ++face) {
                    bool ghostNeighbor = meshReader.getElements()[cell].neighbors[face] == meshReader.getElements().size();
                    const auto& neighbor = [&]() {
                        if (ghostNeighbor) {
                            // ghost layer
                            return container.ghostClusterBackmap.layerPosition(meshReader.getBoundaryElements()[cell].get().boundaryElements[face]);
                        }
                        else {
                            // copy/interior layer
                            return container.clusterBackmap.layerPosition(meshReader.getElements()[cell].neighbors[face]);
                        }
                    }(); // TODO: change here
                    secondaryCellInformation[index].neighborRank[face] = meshReader.getElements()[cell].neighborRanks[face];
                    secondaryCellInformation[index].faceNeighborIds[face] = neighbor.second;
                    cellInformation[index].neighborConfigIds[face] = neighbor.first;
                    cellInformation[index].faceTypes[face] = static_cast<FaceType>(meshReader.getElements()[cell].boundaries[face]);
                    cellInformation[index].faceRelations[face][0] = meshReader.getElements()[cell].neighborSides[face];
                    cellInformation[index].faceRelations[face][1] = meshReader.getElements()[cell].sideOrientations[face];
                    container.cluster.visitIdx(neighbor.first, [&](auto&& ltsview2) {
                        const auto& materialData2 = ltsview2.var(ltsview2.lts.materialData);
                        material[index].neighbor[face] = &materialData2[neighbor.second];
                    });
                }
                ++index;
            }
        }
        else if (layerview.icg == Ghost) {
            std::size_t index = 0;

#ifdef _OPENMP
// #pragma omp parallel for schedule(static)
#endif
            for (const auto& cell : layout.cells(layerview.icg)) {
                new (materialData[index]) MaterialT(materialsDBGhost[cell]);
                secondaryCellInformation[index].meshId = cell;
                secondaryCellInformation[index].configId = layerview.config;
                secondaryCellInformation[index].clusterId = layerview.cluster;

                const auto& boundaryElement = meshReader.getBoundaryElements()[cell].get();
                secondaryCellInformation[index].rank = boundaryElement.neighborRank;
                auto neighbor = container.clusterBackmap.layerPosition(boundaryElement.localElement);
                secondaryCellInformation[index].faceNeighborIds[boundaryElement.localSide] = boundaryElement.localElement;

                int face = boundaryElement.neighborSide;

                secondaryCellInformation[index].neighborRank[face] = MPI::mpi.rank();
                cellInformation[index].neighborConfigIds[face] = configs[boundaryElement.group].configId;
                cellInformation[index].faceTypes[face] = static_cast<FaceType>(meshReader.getElements()[boundaryElement.localElement].boundaries[boundaryElement.localSide]);
                cellInformation[index].faceRelations[face][0] = boundaryElement.localSide;
                cellInformation[index].faceRelations[face][1] = meshReader.getElements()[boundaryElement.localElement].sideOrientations[boundaryElement.localSide];
                ++index;
            }
        }
    });

    // pass 3: LTS setup
    logInfo(rank) << "Setting up LTS configuration...";
    internal::handleLtsSetup(container);

    // pass 4: correct LTS setup, again. Do bucket setup, determine communication datastructures
    logInfo(rank) << "Setting up data exchange and face displacements (buckets)...";
    internal::bucketsAndCommunication(container);

    // TODO init DR and boundary
    logInfo(rank) << "Initializing Dynamic Rupture tree and data structures...";
    allocateDynRupTree(container);

    // allocate boundary tree (effectively behaves like a bucket...)
    logInfo(rank) << "Initializing Boundary tree...";
    allocateBoundaryTree(container);

    logInfo(rank) << "Allocate scratch pads on GPUs...";
    internal::allocateScratchpads(container);

    logInfo(rank) << "Record execution paths on GPUs...";
    internal::allocateScratchpads(container);

    // pass 5: setup cell data
    logInfo(rank) << "Initializing Cell data...";
    initializeCellMatrices(meshReader, container);

    logInfo(rank) << "Setting up GPU-specific data...";
    // TODO:

    // pass 6: time clusters
    logInfo(rank) << "Creating time clusters...";
    container.cluster.visitLayers([&](auto&& layerview) {
        // skip empty time clusters
        if (layerview.layer.getNumberOfCells() > 0) {
            // TODO: add time cluster
            timeManager.addCluster()
        }
    });

    logInfo(rank) << "Connecting time clusters...";
    container.cluster.visitLayers([&](auto&& layerview) {
        // skip empty time clusters
        if (layerview.layer.getNumberOfCells() > 0) {
            // TODO: add time cluster
            timeManager.addCluster()
        }
    });

    logInfo(rank) << "Memory setup done.";
}

} // namespace seissol::initializer
