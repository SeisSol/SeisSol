#include "InputParameters.hpp"

#include "InputAux.hpp"
#include "utils/logger.h"
#include "utils/stringutils.h"

#include "Geometry/MeshReader.h"
#include "SourceTerm/Manager.h"
#include "Checkpoint/Backend.h"
#include "time_stepping/LtsWeights/WeightsFactory.h"
#include "Parallel/MPI.h"

#include <yaml-cpp/yaml.h>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>

using namespace seissol::initializer::parameters;

namespace {

// converts a string to lower case, and trims it.
static void sanitize(std::string& input) {
  utils::StringUtils::trim(input);
  utils::StringUtils::toLower(input);
}

// A small helper class which reads a YAML node dictionary. It keeps track of all items that have
// been read and reports all values which are not used or not used anymore.
// TODO(David): maybe make the reader more tree-like (i.e. keep a central set on which nodes have
// been visited), and output all non-understood values at the end and not between sections
class ParameterReader {
  public:
  ParameterReader(const YAML::Node& node, bool empty) : node(node), empty(empty) {}

  template <typename T>
  T readWithDefault(const std::string& field, const T& defaultValue) {
    T value = defaultValue;
    if (hasField(field)) {
      value = readUnsafe<T>(field);
    } else {
      logDebug(seissol::MPI::mpi.rank())
          << "The field" << field << "was not specified, using fallback.";
    }
    return value;
  }

  // TODO(David): long-term (if we don't switch to another format first), merge readWithDefaultEnum
  // with readWithDefaultStringEnum, i.e. allow both numerical and textual values for an enum (can
  // we maybe auto-generate a parser from an enum definition?)
  template <typename T>
  T readWithDefaultEnum(const std::string& field,
                        const T& defaultValue,
                        const std::unordered_set<T>& validValues) {
    int value = readWithDefault(field, static_cast<int>(defaultValue));
    if (validValues.find(static_cast<T>(value)) == validValues.end()) {
      logError() << "The field" << field << "had an invalid enum value:" << value;
    }
    return static_cast<T>(value);
  }

  template <typename T>
  T readWithDefaultStringEnum(const std::string& field,
                              const std::string& defaultValue,
                              const std::unordered_map<std::string, T>& validValues) {
    std::string value = readWithDefault(field, defaultValue); // TODO(David): sanitize string
    sanitize(value);
    if (validValues.find(value) == validValues.end()) {
      logError() << "The field" << field << "had an invalid enum value:" << value;
    }
    return validValues.at(value);
  }

  template <typename T>
  T readOrFail(const std::string& field, const std::string& failMessage) {
    if (hasField(field)) {
      return readUnsafe<T>(field);
    } else {
      logError() << "The field" << field << "was not found, but it is required.";
      return T(); // unreachable. TODO(David): use compiler hint instead
    }
  }

  void warnDeprecatedSingle(const std::string& field) {
    if (hasField(field)) {
      visited.emplace(field);
      logInfo(seissol::MPI::mpi.rank())
          << "The field" << field
          << "is no longer in use. You may safely remove it from your parameters file.";
    }
  }

  void warnDeprecated(const std::vector<std::string>& fields) {
    for (const auto& field : fields) {
      warnDeprecatedSingle(field);
    }
  }

  void warnUnknown() {
    for (const auto& subnodes : node) {
      auto field = subnodes.first.as<std::string>();
      if (visited.find(field) == visited.end()) {
        logWarning(seissol::MPI::mpi.rank()) << "The field" << field << "is not known to SeisSol.";
      }
    }
  }

  template <typename... Args>
  void markUnused(const Args&... argFields) {
    for (const auto& field : {argFields...}) {
      logDebug(seissol::MPI::mpi.rank()) << "The field" << field << "is ignored (if it is found).";
      visited.emplace(field);
    }
  }

  ParameterReader readSubNode(const std::string& subnodeName) {
    visited.emplace(subnodeName);
    logDebug(seissol::MPI::mpi.rank()) << "Entering section" << subnodeName;
    if (hasField(subnodeName)) {
      return ParameterReader(node[subnodeName], false);
    } else {
      logDebug(seissol::MPI::mpi.rank())
          << "Section" << subnodeName
          << "not found in the given parameter file. Using an empty reader.";
      return ParameterReader(node[subnodeName], true);
    }
  }

  bool hasField(const std::string& field) { return !empty && node[field]; }

  private:
  template <typename T>
  T readUnsafe(const std::string& field) {
    visited.emplace(field);
    logDebug(seissol::MPI::mpi.rank()) << "The field" << field << "was read.";
    try {
      // booleans are stored as integers
      if constexpr (std::is_same<T, bool>::value) {
        return node[field].as<int>() > 0;
      } else {
        return node[field].as<T>();
      }
    } catch (std::exception& e) {
      logError() << "Error while reading field" << field << ":" << e.what();
      return T(); // unreachable. TODO(David): use compiler hint instead
    }
  }

  bool empty;
  YAML::Node node; // apparently the YAML nodes use a reference semantic. Hence, we do it like this.
  std::unordered_set<std::string> visited;
};

static void readModel(ParameterReader& baseReader, SeisSolParameters& seissolParams) {
  auto reader = baseReader.readSubNode("equations");

  seissolParams.model.materialFileName =
      reader.readOrFail<std::string>("materialfilename", "No material file given.");
  seissolParams.model.boundaryFileName =
      reader.readWithDefault("boundaryfilename", std::string(""));
  seissolParams.model.hasBoundaryFile = seissolParams.model.boundaryFileName != "";

  seissolParams.model.gravitationalAcceleration =
      reader.readWithDefault("gravitationalacceleration", 9.81);

  seissolParams.model.plasticity = reader.readWithDefault("plasticity", false);
  seissolParams.model.tv = reader.readWithDefault("tv", 0.1);
  seissolParams.model.useCellHomogenizedMaterial =
      reader.readWithDefault("usecellhomogenizedmaterial", true);
  seissolParams.itmParameters.ITMToggle = reader.readWithDefault("itmenable", bool(0));
  if (seissolParams.itmParameters.ITMToggle) {
    seissolParams.itmParameters.ITMStartingTime = reader.readWithDefault("itmstartingtime", 0.0);
    seissolParams.itmParameters.ITMTime = reader.readWithDefault("itmtime", 0.0);
    seissolParams.itmParameters.ITMVelocityScalingFactor =
        reader.readWithDefault("itmvelocityscalingfactor", 1.0);
    seissolParams.itmParameters.reflectionType =
        reader.readWithDefaultEnum("itmreflectiontype",
                                   ReflectionType::bothwaves,
                                   {ReflectionType::bothwaves,
                                    ReflectionType::bothwaves_velocity,
                                    ReflectionType::pwave,
                                    ReflectionType::swave});
    if (seissolParams.itmParameters.ITMTime <= 0.0) {
      logError() << "ITM Time is not positive. It should be positive!";
    }
    if (seissolParams.itmParameters.ITMVelocityScalingFactor < 0.0) {
      logError() << "ITM Velocity Scaling Factor is less than zero. It should be positive!";
    }
    if (seissolParams.itmParameters.ITMStartingTime < 0.0) {
      logError() << "ITM Starting Time can not be less than zero";
    }
  } else {
    reader.markUnused(
        "itmstartingtime", "itmtime", "itmvelocityscalingfactor", "itmreflectiontype");
  }

  if (isModelViscoelastic()) {
    seissolParams.model.freqCentral = reader.readOrFail<double>(
        "freqcentral", "equations.freqcentral is needed for the attenuation fitting.");
    seissolParams.model.freqRatio = reader.readOrFail<double>(
        "freqratio", "equations.freqratio is needed for the attenuation fitting.");

    if (seissolParams.model.freqRatio <= 0) {
      logError()
          << "The freqratio parameter must be positive---but that is currently not the case.";
    }
  } else {
    reader.markUnused("freqcentral");
    reader.markUnused("freqratio");
  }

  reader.warnDeprecated({"adjoint", "adjfilename", "anisotropy"});
  reader.warnUnknown();
}

static void readCubeGenerator(ParameterReader& baseReader, SeisSolParameters& seissolParams) {
  auto reader = baseReader.readSubNode("cubegenerator");
  seissolParams.cubeGenerator.cubeMinX = reader.readWithDefault("cubeminx", 6);
  seissolParams.cubeGenerator.cubeMaxX = reader.readWithDefault("cubemaxx", 6);
  seissolParams.cubeGenerator.cubeMinY = reader.readWithDefault("cubeminy", 6);
  seissolParams.cubeGenerator.cubeMaxY = reader.readWithDefault("cubemaxy", 6);
  seissolParams.cubeGenerator.cubeMinZ = reader.readWithDefault("cubeminz", 6);
  seissolParams.cubeGenerator.cubeMaxZ = reader.readWithDefault("cubemaxz", 6);

  seissolParams.cubeGenerator.cubeX = reader.readWithDefault("cubex", 2);
  seissolParams.cubeGenerator.cubeY = reader.readWithDefault("cubey", 2);
  seissolParams.cubeGenerator.cubeZ = reader.readWithDefault("cubez", 2);

  // only x dimension has its number of partitions set to number of MPI processes
  seissolParams.cubeGenerator.cubePx = seissol::MPI::mpi.size();
  seissolParams.cubeGenerator.cubePy = 1;
  seissolParams.cubeGenerator.cubePz = 1;

  seissolParams.cubeGenerator.cubeS = reader.readWithDefault("cubes", 100);
  seissolParams.cubeGenerator.cubeSx =
      reader.readWithDefault("cubesx", seissolParams.cubeGenerator.cubeS);
  seissolParams.cubeGenerator.cubeSy =
      reader.readWithDefault("cubesy", seissolParams.cubeGenerator.cubeS);
  seissolParams.cubeGenerator.cubeSz =
      reader.readWithDefault("cubesz", seissolParams.cubeGenerator.cubeS);

  seissolParams.cubeGenerator.cubeTx = reader.readWithDefault("cubetx", 0.0);
  seissolParams.cubeGenerator.cubeTy = reader.readWithDefault("cubety", 0.0);
  seissolParams.cubeGenerator.cubeTz = reader.readWithDefault("cubetz", 0.0);
}

static void readMesh(ParameterReader& baseReader, SeisSolParameters& seissolParams) {
  auto reader = baseReader.readSubNode("meshnml");

  seissolParams.mesh.meshFileName =
      reader.readOrFail<std::string>("meshfile", "No mesh file given.");
  seissolParams.mesh.partitioningLib =
      reader.readWithDefault("partitioninglib", std::string("Default"));
  seissolParams.mesh.meshFormat = reader.readWithDefaultStringEnum<seissol::geometry::MeshFormat>(
      "meshgenerator",
      "puml",
      {{"netcdf", seissol::geometry::MeshFormat::Netcdf},
       {"puml", seissol::geometry::MeshFormat::PUML},
       {"cubegenerator", seissol::geometry::MeshFormat::CubeGenerator}});

  if (seissolParams.mesh.meshFormat == seissol::geometry::MeshFormat::CubeGenerator) {
    readCubeGenerator(baseReader, seissolParams);
  } else {
    reader.markUnused("cubegenerator",
                      "cubeminx",
                      "cubemaxx",
                      "cubeminy",
                      "cubemaxy",
                      "cubeminz",
                      "cubemaxz",
                      "cubex",
                      "cubey",
                      "cubez",
                      "cubes",
                      "cubesx",
                      "cubesy",
                      "cubesz",
                      "cubetx",
                      "cubety",
                      "cubetz");
  }

  seissolParams.mesh.displacement = seissol::initializers::convertStringToArray<double, 3>(
      reader.readWithDefault("displacement", std::string("0.0 0.0 0.0")));
  auto scalingX = seissol::initializers::convertStringToArray<double, 3>(
      reader.readWithDefault("scalingmatrixx", std::string("1.0 0.0 0.0")));
  auto scalingY = seissol::initializers::convertStringToArray<double, 3>(
      reader.readWithDefault("scalingmatrixy", std::string("0.0 1.0 0.0")));
  auto scalingZ = seissol::initializers::convertStringToArray<double, 3>(
      reader.readWithDefault("scalingmatrixz", std::string("0.0 0.0 1.0")));
  seissolParams.mesh.scaling = {scalingX, scalingY, scalingZ};

  seissolParams.timeStepping.vertexWeight.weightElement =
      reader.readWithDefault("vertexweightelement", 100);
  seissolParams.timeStepping.vertexWeight.weightDynamicRupture =
      reader.readWithDefault("vertexweightdynamicrupture", 100);
  seissolParams.timeStepping.vertexWeight.weightFreeSurfaceWithGravity =
      reader.readWithDefault("vertexweightfreesurfacewithgravity", 100);

  seissolParams.mesh.showEdgeCutStatistics = reader.readWithDefault("showedgecutstatistics", false);

  reader.warnDeprecated({"periodic", "periodic_direction"});
  reader.warnUnknown();
}

static void readTimeStepping(ParameterReader& baseReader, SeisSolParameters& seissolParams) {
  auto reader = baseReader.readSubNode("discretization");

  seissolParams.timeStepping.cfl = reader.readWithDefault("cfl", 0.5);
  seissolParams.timeStepping.lts.rate = reader.readWithDefault("clusteredlts", 2u);
  seissolParams.timeStepping.lts.weighttype = reader.readWithDefaultEnum(
      "ltsweighttypeid",
      seissol::initializers::time_stepping::LtsWeightsTypes::ExponentialWeights,
      {
          seissol::initializers::time_stepping::LtsWeightsTypes::ExponentialWeights,
          seissol::initializers::time_stepping::LtsWeightsTypes::ExponentialBalancedWeights,
          seissol::initializers::time_stepping::LtsWeightsTypes::EncodedBalancedWeights,
      });

  if (isModelViscoelastic()) {
    // NOTE: we are using a half-initialized struct here... (i.e. be careful)
    double maxTimestepWidthDefault =
        0.25 / (seissolParams.model.freqCentral * std::sqrt(seissolParams.model.freqRatio));
    if (reader.hasField("fixtimestep")) {
      seissolParams.timeStepping.maxTimestepWidth =
          reader.readWithDefault("fixtimestep", maxTimestepWidthDefault);
      if (seissolParams.timeStepping.maxTimestepWidth > maxTimestepWidthDefault) {
        logWarning(seissol::MPI::mpi.rank())
            << "The given maximum timestep width (fixtimestep) is set to"
            << seissolParams.timeStepping.maxTimestepWidth
            << "which is larger than the recommended value of" << maxTimestepWidthDefault
            << " for visco-elastic material (as specified in the documentation). Please be aware "
               "that a too large maximum timestep width may cause the solution to become unstable.";
      } else {
        logInfo(seissol::MPI::mpi.rank())
            << "Maximum timestep width (fixtimestep) given as"
            << seissolParams.timeStepping.maxTimestepWidth << "(less or equal to reference timestep"
            << maxTimestepWidthDefault << ")";
      }
    } else {
      seissolParams.timeStepping.maxTimestepWidth = maxTimestepWidthDefault;
      logInfo(seissol::MPI::mpi.rank())
          << "Setting maximum timestep width to" << maxTimestepWidthDefault
          << " for visco-elastic material (as specified in the documentation).";
    }
  } else {
    seissolParams.timeStepping.maxTimestepWidth = reader.readWithDefault("fixtimestep", 5000.0);
  }

  // TODO(David): integrate LTS parameters here
  reader.markUnused("ltswigglefactormin",
                    "ltswigglefactorstepsize",
                    "ltswigglefactorenforcemaximumdifference",
                    "ltsmaxnumberofclusters",
                    "ltsautomergeclusters",
                    "ltsallowedrelativeperformancelossautomerge",
                    "ltsautomergecostbaseline");

  reader.warnDeprecated({"ckmethod",
                         "dgfineout1d",
                         "fluxmethod",
                         "iterationcriterion",
                         "npoly",
                         "npolyrec",
                         "limitersecurityfactor",
                         "order",
                         "material",
                         "npolymap"});
  reader.warnUnknown();
}

static void readInitialization(ParameterReader& baseReader, SeisSolParameters& seissolParams) {
  auto reader = baseReader.readSubNode("inicondition");

  seissolParams.initialization.type = reader.readWithDefaultStringEnum<InitializationType>(
      "cictype",
      "zero",
      {
          {"zero", InitializationType::Zero},
          {"planarwave", InitializationType::Planarwave},
          {"superimposedplanarwave", InitializationType::SuperimposedPlanarwave},
          {"travelling", InitializationType::Travelling},
          {"acoustictravellingwithitm", InitializationType::AcousticTravellingwithITM},
          {"scholte", InitializationType::Scholte},
          {"snell", InitializationType::Snell},
          {"ocean_0", InitializationType::Ocean0},
          {"ocean_1", InitializationType::Ocean1},
          {"ocean_2", InitializationType::Ocean2},
          {"pressureinjection", InitializationType::PressureInjection},
      });
  const auto originString = reader.readWithDefault("origin", std::string("0.0 0.0 0.0"));
  seissolParams.initialization.origin =
      seissol::initializers::convertStringToArray<double, 3>(originString);
  const auto kVecString = reader.readWithDefault("kvec", std::string("0.0 0.0 0.0"));
  seissolParams.initialization.kVec =
      seissol::initializers::convertStringToArray<double, 3>(kVecString);
  std::string defaultAmpFieldString;
  for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
    defaultAmpFieldString += " 0.0";
  }
  const auto ampFieldString = reader.readWithDefault("ampfield", defaultAmpFieldString);
  seissolParams.initialization.ampField =
      seissol::initializers::convertStringToArray<double, NUMBER_OF_QUANTITIES>(ampFieldString);
  seissolParams.initialization.magnitude = reader.readWithDefault("magnitude", 0.0);
  seissolParams.initialization.width =
      reader.readWithDefault("width", std::numeric_limits<double>::infinity());
  seissolParams.initialization.k = reader.readWithDefault("k", 0.0);

  reader.warnUnknown();
}

static void readOutput(ParameterReader& baseReader, SeisSolParameters& seissolParams) {
  auto reader = baseReader.readSubNode("output");

  constexpr double veryLongTime = 1.0e100;

  auto warnIntervalAndDisable =
      [](bool& enabled, double interval, const std::string& valName, const std::string& intName) {
        if (enabled && interval <= 0) {
          auto intPhrase = valName + " = 0";
          logInfo(seissol::MPI::mpi.rank())
              << "In your parameter file, you have specified a non-positive interval for" << intName
              << ". This mechanism is deprecated and may be removed in a future version of "
                 "SeisSol. Consider disabling the whole module by setting"
              << valName << "to 0 instead by adding" << intPhrase
              << "to the \"output\" section of your parameter file instead.";

          // still, replicate the old behavior.
          enabled = false;
        }
      };

  // general params
  seissolParams.output.format = reader.readWithDefaultEnum<OutputFormat>(
      "format", OutputFormat::None, {OutputFormat::None, OutputFormat::Xdmf});
  seissolParams.output.prefix =
      reader.readOrFail<std::string>("outputfile", "Output file prefix not defined.");
  seissolParams.output.xdmfWriterBackend =
      reader.readWithDefaultStringEnum<xdmfwriter::BackendType>(
          "xdmfwriterbackend",
          "posix",
          {
              {"posix", xdmfwriter::BackendType::POSIX},
#ifdef USE_HDF
              {"hdf5", xdmfwriter::BackendType::H5},
#endif
          });

  // checkpointing
  seissolParams.output.checkpointParameters.enabled = reader.readWithDefault("checkpoint", true);
  seissolParams.output.checkpointParameters.backend =
      reader.readWithDefaultStringEnum<seissol::checkpoint::Backend>(
          "checkpointbackend",
          "none",
          {{"none", seissol::checkpoint::Backend::DISABLED},
           {"posix", seissol::checkpoint::Backend::POSIX},
           {"hdf5", seissol::checkpoint::Backend::HDF5},
           {"mpio", seissol::checkpoint::Backend::MPIO},
           {"mpio_async", seissol::checkpoint::Backend::MPIO_ASYNC},
           {"sionlib", seissol::checkpoint::Backend::SIONLIB}});
  seissolParams.output.checkpointParameters.interval =
      reader.readWithDefault("checkpointinterval", 0.0);

  warnIntervalAndDisable(seissolParams.output.checkpointParameters.enabled,
                         seissolParams.output.checkpointParameters.interval,
                         "checkpoint",
                         "checkpointinterval");

  if (seissolParams.output.checkpointParameters.enabled) {
    seissolParams.output.checkpointParameters.fileName =
        reader.readOrFail<std::string>("checkpointfile", "No checkpoint filename given.");
  } else {
    reader.markUnused("checkpointfile");
  }

  // output: wavefield
  // (these variables are usually not prefixed with "wavefield" or the likes)

  // bounds
  auto bounds = seissol::initializers::convertStringToArray<double, 6>(
      reader.readWithDefault("outputregionbounds", std::string("0.0 0.0 0.0 0.0 0.0 0.0")));
  seissolParams.output.waveFieldParameters.bounds.boundsX.lower = bounds[0];
  seissolParams.output.waveFieldParameters.bounds.boundsX.upper = bounds[1];
  seissolParams.output.waveFieldParameters.bounds.boundsY.lower = bounds[2];
  seissolParams.output.waveFieldParameters.bounds.boundsY.upper = bounds[3];
  seissolParams.output.waveFieldParameters.bounds.boundsZ.lower = bounds[4];
  seissolParams.output.waveFieldParameters.bounds.boundsZ.upper = bounds[5];
  seissolParams.output.waveFieldParameters.bounds.enabled =
      !(bounds[0] == 0 && bounds[1] == 0 && bounds[2] == 0 && bounds[3] == 0 && bounds[4] == 0 &&
        bounds[5] == 0);

  seissolParams.output.waveFieldParameters.enabled =
      reader.readWithDefault("wavefieldoutput", true);
  seissolParams.output.waveFieldParameters.interval =
      reader.readWithDefault("timeinterval", veryLongTime);
  seissolParams.output.waveFieldParameters.refinement =
      reader.readWithDefaultEnum<OutputRefinement>("refinement",
                                                   OutputRefinement::NoRefine,
                                                   {OutputRefinement::NoRefine,
                                                    OutputRefinement::Refine4,
                                                    OutputRefinement::Refine8,
                                                    OutputRefinement::Refine32});

  warnIntervalAndDisable(seissolParams.output.waveFieldParameters.enabled,
                         seissolParams.output.waveFieldParameters.interval,
                         "wavefieldoutput",
                         "timeinterval");

  if (seissolParams.output.waveFieldParameters.enabled &&
      seissolParams.output.format == OutputFormat::None) {
    logInfo(seissol::MPI::mpi.rank())
        << "Disabling the wavefield output by setting \"outputformat = 10\" is deprecated "
           "and may be removed in a future version of SeisSol. Consider using the parameter "
           "\"wavefieldoutput\" instead. To disable wavefield output, add \"wavefieldoutput "
           "= 0\" to the \"output\" section of your parameters file.";

    seissolParams.output.waveFieldParameters.enabled = false;
  }

  auto groupsVector = reader.readWithDefault("outputgroups", std::vector<int>());
  seissolParams.output.waveFieldParameters.groups =
      std::unordered_set<int>(groupsVector.begin(), groupsVector.end());

  // output mask
  auto iOutputMask = reader.readOrFail<std::string>("ioutputmask", "No output mask given.");
  seissol::initializers::convertStringToMask(iOutputMask,
                                             seissolParams.output.waveFieldParameters.outputMask);

  auto iPlasticityMask = reader.readWithDefault("iplasticitymask", std::string("0 0 0 0 0 0 1"));
  seissol::initializers::convertStringToMask(
      iPlasticityMask, seissolParams.output.waveFieldParameters.plasticityMask);

  auto integrationMask =
      reader.readWithDefault("integrationmask", std::string("0 0 0 0 0 0 0 0 0"));
  seissol::initializers::convertStringToMask(
      integrationMask, seissolParams.output.waveFieldParameters.integrationMask);

  // output: surface
  seissolParams.output.freeSurfaceParameters.enabled =
      reader.readWithDefault("surfaceoutput", false);
  seissolParams.output.freeSurfaceParameters.interval =
      reader.readWithDefault("surfaceoutputinterval", veryLongTime);
  seissolParams.output.freeSurfaceParameters.refinement =
      reader.readWithDefault("surfaceoutputrefinement", 0u);

  warnIntervalAndDisable(seissolParams.output.freeSurfaceParameters.enabled,
                         seissolParams.output.freeSurfaceParameters.interval,
                         "surfaceoutput",
                         "surfaceoutputinterval");

  // output: energy
  seissolParams.output.energyParameters.enabled = reader.readWithDefault("energyoutput", false);
  seissolParams.output.energyParameters.interval =
      reader.readWithDefault("energyoutputinterval", veryLongTime);
  seissolParams.output.energyParameters.terminalOutput =
      reader.readWithDefault("energyterminaloutput", false);
  seissolParams.output.energyParameters.computeVolumeEnergiesEveryOutput =
      reader.readWithDefault("computevolumeenergieseveryoutput", 1);

  warnIntervalAndDisable(seissolParams.output.energyParameters.enabled,
                         seissolParams.output.energyParameters.interval,
                         "energyoutput",
                         "energyoutputinterval");

  // output: refinement
  seissolParams.output.receiverParameters.enabled = reader.readWithDefault("receiveroutput", true);
  seissolParams.output.receiverParameters.interval =
      reader.readWithDefault("receiveroutputinterval", veryLongTime);
  seissolParams.output.receiverParameters.computeRotation =
      reader.readWithDefault("receivercomputerotation", false);
  seissolParams.output.receiverParameters.fileName =
      reader.readWithDefault("rfilename", std::string(""));
  seissolParams.output.receiverParameters.samplingInterval = reader.readWithDefault("pickdt", 0.0);

  warnIntervalAndDisable(seissolParams.output.receiverParameters.enabled,
                         seissolParams.output.receiverParameters.interval,
                         "receiveroutput",
                         "receiveroutputinterval");

  // output: loop statistics
  seissolParams.output.loopStatisticsNetcdfOutput =
      reader.readWithDefault("loopstatisticsnetcdfoutput", false);

  reader.warnDeprecated({"rotation",
                         "interval",
                         "nrecordpoints",
                         "printintervalcriterion",
                         "pickdttype",
                         "ioutputmaskmaterial",
                         "faultoutputflag"});
  reader.warnUnknown();
}

static void readAbortCriteria(ParameterReader& baseReader, SeisSolParameters& seissolParams) {
  auto reader = baseReader.readSubNode("abortcriteria");

  seissolParams.end.endTime = reader.readWithDefault("endtime", 15.0);

  reader.warnDeprecated(
      {"maxiterations", "maxtolerance", "maxtolcriterion", "walltime_h", "delay_h"});
  reader.warnUnknown();
}

static void readSource(ParameterReader& baseReader, SeisSolParameters& seissolParams) {
  auto reader = baseReader.readSubNode("sourcetype");

  seissolParams.source.type =
      reader.readWithDefaultEnum("type",
                                 seissol::sourceterm::SourceType::None,
                                 {seissol::sourceterm::SourceType::None,
                                  seissol::sourceterm::SourceType::FsrmSource,
                                  seissol::sourceterm::SourceType::NrfSource});
  if (seissolParams.source.type != seissol::sourceterm::SourceType::None) {
    seissolParams.source.fileName =
        reader.readOrFail<std::string>("filename", "No source file specified.");
  } else {
    reader.markUnused("filename");
  }

  reader.warnDeprecated({"rtype", "ndirac", "npulsesource", "nricker"});
  reader.warnUnknown();
}

} // namespace

void SeisSolParameters::readParameters(const YAML::Node& baseNode) {
  logInfo(seissol::MPI::mpi.rank()) << "Reading SeisSol parameter file...";

  ParameterReader baseReader(baseNode, false);

  readModel(baseReader, *this);
  readMesh(baseReader, *this);
  readTimeStepping(baseReader, *this);
  readInitialization(baseReader, *this);
  readOutput(baseReader, *this);
  readSource(baseReader, *this);
  readAbortCriteria(baseReader, *this);

  // TODO(David): remove once DR parameter reading is integrated here
  baseReader.markUnused("dynamicrupture", "elementwise", "pickpoint");

  baseReader.warnDeprecated({"boundaries",
                             "rffile",
                             "inflowbound",
                             "inflowboundpwfile",
                             "inflowbounduin",
                             "source110",
                             "source15",
                             "source1618",
                             "source17",
                             "source19",
                             "spongelayer",
                             "sponges",
                             "analysis",
                             "analysisfields",
                             "debugging"});
  baseReader.warnUnknown();

  logInfo(seissol::MPI::mpi.rank()) << "SeisSol parameter file read successfully.";

  auto printYesNo = [](bool yesno) { return yesno ? "yes" : "no"; };

  logInfo(seissol::MPI::mpi.rank()) << "Model information:";
  logInfo(seissol::MPI::mpi.rank()) << "Elastic model:" << printYesNo(isModelElastic());
  logInfo(seissol::MPI::mpi.rank()) << "Viscoelastic model:" << printYesNo(isModelViscoelastic());
  logInfo(seissol::MPI::mpi.rank()) << "Anelastic model:" << printYesNo(isModelAnelastic());
  logInfo(seissol::MPI::mpi.rank()) << "Poroelastic model:" << printYesNo(isModelPoroelastic());
  logInfo(seissol::MPI::mpi.rank()) << "Anisotropic model:" << printYesNo(isModelAnisotropic());
  logInfo(seissol::MPI::mpi.rank()) << "Plasticity:" << printYesNo(this->model.plasticity);
}
