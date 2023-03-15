
#ifndef INPUT_PARAMETERS_HPP_
#define INPUT_PARAMETERS_HPP_

#include <cstdint>
#include <string>
#include <array>
#include <yaml-cpp/yaml.h>

#include "Geometry/MeshReader.h"
#include "SourceTerm/Manager.h"

namespace seissol::initializer::parameters {
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
        double freqCentral;
        double freqRatio;
        std::string materialFileName = "";
        std::string boundaryFileName = "";

        void readPar(const YAML::Node& baseNode);
    };

    struct InitializationParameters {
        std::string type;
        std::array<double, 3> origin;
        std::array<double, 3> kVec;
        std::array<double, NUMBER_OF_QUANTITIES> ampField;
    };

    struct DynamicRuptureParameters {
        bool hasFault;
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

    struct OutputReceiverParameters {
    };

    struct OutputSurfaceParameters {

    };

    struct OutputEnergyParameters {

    };

    struct OutputCheckpointParameters {
        double interval;
        std::string fileName;
        std::string backend;
    };

    struct OutputInterval {
        double lower;
        double higher;
    };

    struct OutputBounds {
        OutputInterval boundsX, boundsY, boundsZ;
    };

    struct OutputParameters {
        std::string prefix = "data";
        OutputFormat format = OutputFormat::None;
        OutputRefinement refinement = OutputRefinement::NoRefine;
        OutputBounds outputBounds;
        std::vector<bool> outputMask;
        std::vector<bool> plasticityMask;
        std::vector<bool> integrationMask;
        std::string xdmfWriterBackend;
    };

    struct LtsParameters {
        unsigned rate = 2;
        int weighttype = 1; // TODO: type
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
