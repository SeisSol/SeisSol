#ifndef FACECLUSTER_H_
#define FACECLUSTER_H_

#include "AbstractTimeCluster.h"
#include "Kernels/DynamicRupture.h"
#include <DynamicRupture/FrictionLaws/FrictionSolver.h>
#include <DynamicRupture/Output/OutputManager.hpp>
#include <Initializer/DynamicRupture.h>
#include <Initializer/tree/Layer.hpp>
#include <Monitoring/ActorStateStatistics.h>
#include <Parallel/Runtime/Stream.hpp>
#include <Solver/time_stepping/ActorState.h>

namespace seissol::time_stepping {
    class DynamicRuptureCluster : public FaceCluster
    {
    public:
    DynamicRuptureCluster(double maxTimeStepSize,
                                                 long timeStepRate,
                                                 unsigned profilingId,
        seissol::initializer::Layer* layer,
        seissol::initializer::DynamicRupture* descr,
        CompoundGlobalData globalData,
        seissol::dr::friction_law::FrictionSolver* frictionSolver,
        seissol::dr::friction_law::FrictionSolver* frictionSolverDevice,
        seissol::dr::output::OutputManager* faultOutputManager,
        seissol::SeisSol& seissolInstance,
        LoopStatistics* loopStatistics,
        ActorStateStatistics* actorStateStatistics)
        : FaceCluster(
            maxTimeStepSize, timeStepRate,
#ifdef ACL_DEVICE
      layer->getNumberOfCells() >= deviceHostSwitch() ? Executor::Device : Executor::Host
#else
      Executor::Host
#endif
        ), profilingId(profilingId),
        layer(layer), descr(descr), frictionSolver(frictionSolver), frictionSolverDevice(frictionSolverDevice), faultOutputManager(faultOutputManager),
        globalDataOnHost(globalData.onHost), globalDataOnDevice(globalData.onDevice), loopStatistics(loopStatistics), actorStateStatistics(actorStateStatistics),
        seissolInstance(seissolInstance) {
        dynamicRuptureKernel.setGlobalData(globalData);
        // TODO(David): activate here already
        // frictionSolver->copyLtsTreeToLocal(*layer, descr, 0);
        // frictionSolverDevice->copyLtsTreeToLocal(*layer, descr, 0);
        regionComputeDynamicRupture = loopStatistics->getRegion("computeDynamicRupture");
    }

    void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream);
    void computeFlops();

    void setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager) {
        // TODO(David): remove
        this->faultOutputManager = faultOutputManager;
    }
    std::string description() const override {
        return "face-cluster";
    }
    ~DynamicRuptureCluster() override = default;
    protected:
    seissol::SeisSol& seissolInstance;

    void computeDynamicRupture(  );
    void handleDynamicRupture(  );
    void computeDynamicRuptureDevice(  );

    void computeDynamicRuptureFlops(long long& nonZeroFlops,
                                    long long& hardwareFlops);
    
    void handleAdvancedComputeTimeMessage(ComputeStep step, const NeighborCluster& neighborCluster) override {}
    void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override {}

    seissol::initializer::Layer* layer;
    seissol::initializer::DynamicRupture* descr;
    seissol::dr::friction_law::FrictionSolver* frictionSolver;
    seissol::dr::friction_law::FrictionSolver* frictionSolverDevice;
    seissol::dr::output::OutputManager* faultOutputManager;

    LoopStatistics* loopStatistics;
    ActorStateStatistics* actorStateStatistics;
    unsigned        regionComputeDynamicRupture;

    seissol::kernels::DynamicRupture dynamicRuptureKernel;

    GlobalData *globalDataOnHost{nullptr};
    GlobalData *globalDataOnDevice{nullptr};
#ifdef ACL_DEVICE
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

    bool computePickpoint{false};

    ClusterTimes ct;

    long long flops_nonZero;
    long long flops_hardware;

    unsigned profilingId;

    void start() override {}

    void runCompute(ComputeStep step) override;
    };
} // namespace seissol::time_stepping

#endif
