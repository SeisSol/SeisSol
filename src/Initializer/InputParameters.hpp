
#ifndef INPUT_PARAMETERS_HPP_
#define INPUT_PARAMETERS_HPP_

#include <cstdint>
#include <string>
#include <array>
#include <yaml-cpp/yaml.h>

#include <xdmfwriter/XdmfWriter.h>

#include "Geometry/MeshReader.h"
#include "SourceTerm/Manager.h"
#include "Checkpoint/Backend.h"
#include "time_stepping/LtsWeights/WeightsFactory.h"

namespace seissol::initializer::parameters {
    //constexpr auto NUMBER_OF_QUANTITIES = tensor::Q::Shape[ sizeof(tensor::Q::Shape) / sizeof(tensor::Q::Shape[0]) - 1];

    constexpr bool modelAnelastic() {
        return NUMBER_OF_RELAXATION_MECHANISMS > 0;
    }

    constexpr bool modelPoroelastic() {
#ifdef USE_POROELASTIC
        return true;
#else
        return false;
#endif
    }

    constexpr bool modelAnisotropic() {
#ifdef USE_ANISOTROPIC
        return true;
#else
        return false;
#endif
    }

    struct ModelParameters {
        double gravitationalAcceleration = 9.81;
        double tv = 0.1;
        bool plasticity = false;
        bool useCellHomogenizedMaterial = false;
        double freqCentral = 0.0;
        double freqRatio = 0.0;
        std::string materialFileName = "";
        std::string boundaryFileName = "";
    };

    enum class InitializationType {
        Zero,
        Planarwave,
        SuperimposedPlanarwave,
        Travelling,
        Scholte,
        Snell,
        Ocean0,
        Ocean1,
        Ocean2
    };

    struct InitializationParameters {
        InitializationType type;
        std::array<double, 3> origin;
        std::array<double, 3> kVec;
        std::array<double, NUMBER_OF_QUANTITIES> ampField;
    };

    struct DynamicRuptureParameters {
        bool hasFault;
        // TODO: port rest of the DR parameters here?
    };

    enum class OutputFormat {
        None = 10,
        Xdmf = 6
    };

    enum class OutputRefinement {
        NoRefine = 0,
        Refine4 = 1,
        Refine8 = 2,
        Refine32 = 3
    };

    struct VertexWeightParameters{
        int weightElement;
        int weightDynamicRupture;
        int weightFreeSurfaceWithGravity;
    };

    struct MeshParameters {
        std::string meshFileName;
        MeshFormat meshFormat;
        std::array<double, 3> displacement;
        std::array<std::array<double, 3>, 3> scaling;
    };

    struct OutputInterval {
        double lower;
        double upper;

        bool contains(double value) const {
            return value >= lower && value <= upper;
        }
    };

    struct OutputBounds {
        bool enabled = false;
        OutputInterval boundsX, boundsY, boundsZ;

        bool contains(double x, double y, double z) const {
            if (enabled) {
                return boundsX.contains(x) && boundsY.contains(y) && boundsZ.contains(z);
            }
            else {
                return true;
            }
        }
    };

    struct OutputBaseParameters {
        bool enabled = false;
        double interval = 1.0e100;
    };

    struct OutputReceiverParameters : public OutputBaseParameters {
        bool computeRotation = false;
        std::string fileName;
        double samplingInterval = 0.0;
    };

    struct OutputSurfaceParameters : public OutputBaseParameters {
        unsigned refinement;
    };

    struct OutputEnergyParameters : public OutputBaseParameters {
        bool terminalOutput = false;
        bool computeVolumeEnergiesEveryOutput = true;
    };

    struct CheckpointParameters : public OutputBaseParameters {
        std::string fileName;
        seissol::checkpoint::Backend backend;
    };

    struct OutputWaveFieldParameters : public OutputBaseParameters {
        OutputRefinement refinement = OutputRefinement::NoRefine;
        OutputBounds bounds;
        std::array<bool, NUMBER_OF_QUANTITIES> outputMask;
        std::array<bool, 7> plasticityMask;
        std::array<bool, 9> integrationMask;
        std::unordered_set<int> groups;
    };

    struct OutputParameters {
        std::string prefix = "data";
        OutputFormat format = OutputFormat::None;
        xdmfwriter::BackendType xdmfWriterBackend;
        CheckpointParameters checkpointParameters;
        OutputWaveFieldParameters waveFieldParameters;
        OutputReceiverParameters receiverParameters;
        OutputSurfaceParameters freeSurfaceParameters;
        OutputEnergyParameters energyParameters;
        bool faultOutput = false;
        bool loopStatisticsNetcdfOutput = false;
    };

    struct LtsParameters {
        unsigned rate = 2;
        seissol::initializers::time_stepping::LtsWeightsTypes weighttype = seissol::initializers::time_stepping::LtsWeightsTypes::ExponentialWeights;
    };

    struct TimesteppingParameters {
        double cfl = 0.5;
        double maxTimestep = 5000;
        LtsParameters lts;
        VertexWeightParameters vertexWeight;
    };

    struct SourceParameters {
        seissol::sourceterm::SourceType type = seissol::sourceterm::SourceType::None;
        std::string fileName = "";
    };

    struct EndParameters {
        double endTime = 15.0;
        uint64_t maxIterations = 100000000;
    };

    struct SeisSolParameters {
        ModelParameters model;
        DynamicRuptureParameters dynamicRupture;
        MeshParameters mesh;
        InitializationParameters initialization;
        OutputParameters output;
        TimesteppingParameters timestepping;
        SourceParameters source;
        EndParameters end;

        void readPar(const YAML::Node& baseNode);
        void printInfo();
    };
}

#endif
