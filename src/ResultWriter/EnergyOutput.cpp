// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EnergyOutput.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "Equations/Energy.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/PreProcessorMacros.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Modules/Modules.h"
#include "Monitoring/Unit.h"
#include "Numerical/Quadrature.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <ios>
#include <limits>
#include <mpi.h>
#include <optional>
#include <ostream>
#include <string>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

GENERATE_HAS_MEMBER(vInv)
GENERATE_HAS_MEMBER(evalAtQP)

namespace seissol::writer {

namespace {

std::array<real, multisim::NumSimulations>
    computeStaticWork(const real* degreesOfFreedomPlus,
                      const real* degreesOfFreedomMinus,
                      const DRFaceInformation& faceInfo,
                      const DRGodunovData& godunovData,
                      const real slip[seissol::tensor::slipInterpolated::size()],
                      const GlobalData* global) {
  real points[seissol::kernels::NumSpaceQuadraturePoints][2];
  alignas(Alignment) real spaceWeights[seissol::kernels::NumSpaceQuadraturePoints];
  seissol::quadrature::TriangleQuadrature(points, spaceWeights, ConvergenceOrder + 1);

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints krnl;
  krnl.V3mTo2n = global->faceToNodalMatrices;

  alignas(PagesizeStack) real qInterpolatedPlus[tensor::QInterpolatedPlus::size()];
  alignas(PagesizeStack) real qInterpolatedMinus[tensor::QInterpolatedMinus::size()];
  alignas(Alignment) real tractionInterpolated[tensor::tractionInterpolated::size()];
  alignas(Alignment) real qPlus[tensor::Q::size()];
  alignas(Alignment) real qMinus[tensor::Q::size()];

  // needed to counter potential mis-alignment
  std::memcpy(qPlus, degreesOfFreedomPlus, sizeof(qPlus));
  std::memcpy(qMinus, degreesOfFreedomMinus, sizeof(qMinus));

  krnl.QInterpolated = qInterpolatedPlus;
  krnl.Q = qPlus;
  krnl.TinvT = godunovData.dataTinvT;
  krnl._prefetch.QInterpolated = qInterpolatedPlus;
  krnl.execute(faceInfo.plusSide, 0);

  krnl.QInterpolated = qInterpolatedMinus;
  krnl.Q = qMinus;
  krnl.TinvT = godunovData.dataTinvT;
  krnl._prefetch.QInterpolated = qInterpolatedMinus;
  krnl.execute(faceInfo.minusSide, faceInfo.faceRelation);

  dynamicRupture::kernel::computeTractionInterpolated trKrnl;
  trKrnl.tractionPlusMatrix = godunovData.tractionPlusMatrix;
  trKrnl.tractionMinusMatrix = godunovData.tractionMinusMatrix;
  trKrnl.QInterpolatedPlus = qInterpolatedPlus;
  trKrnl.QInterpolatedMinus = qInterpolatedMinus;
  trKrnl.tractionInterpolated = tractionInterpolated;
  trKrnl.execute();

  alignas(Alignment) real staticFrictionalWork[tensor::staticFrictionalWork::size()]{};

  dynamicRupture::kernel::accumulateStaticFrictionalWork feKrnl;
  feKrnl.slipInterpolated = slip;
  feKrnl.tractionInterpolated = tractionInterpolated;
  feKrnl.spaceWeights = spaceWeights;
  feKrnl.staticFrictionalWork = staticFrictionalWork;
  feKrnl.minusSurfaceArea = -0.5 * godunovData.doubledSurfaceArea;
  feKrnl.execute();

  std::array<real, multisim::NumSimulations> frictionalWorkReturn{};
  std::copy_n(staticFrictionalWork, multisim::NumSimulations, frictionalWorkReturn.begin());
  return frictionalWorkReturn;
}

} // namespace

double& EnergiesStorage::gravitationalEnergy(size_t sim) {
  return energies[0 + sim * NumberOfEnergies];
}
double& EnergiesStorage::acousticEnergy(size_t sim) { return energies[1 + sim * NumberOfEnergies]; }
double& EnergiesStorage::acousticKineticEnergy(size_t sim) {
  return energies[2 + sim * NumberOfEnergies];
}
double& EnergiesStorage::elasticEnergy(size_t sim) { return energies[3 + sim * NumberOfEnergies]; }
double& EnergiesStorage::elasticKineticEnergy(size_t sim) {
  return energies[4 + sim * NumberOfEnergies];
}
double& EnergiesStorage::totalFrictionalWork(size_t sim) {
  return energies[5 + sim * NumberOfEnergies];
}
double& EnergiesStorage::staticFrictionalWork(size_t sim) {
  return energies[6 + sim * NumberOfEnergies];
}
double& EnergiesStorage::plasticMoment(size_t sim) { return energies[7 + sim * NumberOfEnergies]; }
double& EnergiesStorage::seismicMoment(size_t sim) { return energies[8 + sim * NumberOfEnergies]; }
double& EnergiesStorage::potency(size_t sim) { return energies[9 + sim * NumberOfEnergies]; }

double& EnergiesStorage::totalMomentumX(size_t sim) {
  return energies[10 + sim * NumberOfEnergies];
}
double& EnergiesStorage::totalMomentumY(size_t sim) {
  return energies[11 + sim * NumberOfEnergies];
}
double& EnergiesStorage::totalMomentumZ(size_t sim) {
  return energies[12 + sim * NumberOfEnergies];
}

void EnergyOutput::init(
    GlobalData* newGlobal,
    const DynamicRupture::Storage& newDynRuptTree,
    const seissol::geometry::MeshReader& newMeshReader,
    const LTS::Storage& newStorage,
    bool newIsPlasticityEnabled,
    const std::string& outputFileNamePrefix,
    const seissol::initializer::parameters::EnergyOutputParameters& parameters) {
  if (parameters.enabled && parameters.interval > 0) {
    isEnabled_ = true;
  } else {
    return;
  }
  const auto rank = Mpi::mpi.rank();
  logInfo() << "Initializing energy output.";

  energyOutputInterval_ = parameters.interval;
  isFileOutputEnabled_ = rank == 0;
  isTerminalOutputEnabled_ = parameters.terminalOutput && (rank == 0);
  terminatorMaxTimePostRupture_ = parameters.terminatorMaxTimePostRupture;
  terminatorMomentRateThreshold_ = parameters.terminatorMomentRateThreshold;
  isCheckAbortCriteraSlipRateEnabled_ = std::isfinite(terminatorMaxTimePostRupture_);
  isCheckAbortCriteraMomentRateEnabled_ = (terminatorMomentRateThreshold_ > 0);
  computeVolumeEnergiesEveryOutput_ = parameters.computeVolumeEnergiesEveryOutput;
  outputFileName_ = outputFileNamePrefix + "-energy.csv";

  global_ = newGlobal;
  drStorage_ = &newDynRuptTree;
  meshReader_ = &newMeshReader;
  ltsStorage_ = &newStorage;

  isPlasticityEnabled_ = newIsPlasticityEnabled;

  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  setSyncInterval(parameters.interval);
}

void EnergyOutput::syncPoint(double time) {
  assert(isEnabled_);
  const auto rank = Mpi::mpi.rank();
  logInfo() << "Writing energy output at time" << time;

  seissolInstance_.dofSync().syncDofs(time);

  computeEnergies();
  reduceEnergies();
  if (isCheckAbortCriteraSlipRateEnabled_) {
    reduceMinTimeSinceSlipRateBelowThreshold();
  }
  if ((rank == 0) && isCheckAbortCriteraMomentRateEnabled_) {
    for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
      const double seismicMomentRate =
          (energiesStorage_.seismicMoment(sim) - seismicMomentPrevious_[sim]) /
          energyOutputInterval_;
      seismicMomentPrevious_[sim] = energiesStorage_.seismicMoment(sim);
      if (time > 0 && seismicMomentRate < terminatorMomentRateThreshold_) {
        minTimeSinceMomentRateBelowThreshold_[sim] += energyOutputInterval_;
      } else {
        minTimeSinceMomentRateBelowThreshold_[sim] = 0.0;
      }
    }
  }
  if (isTerminalOutputEnabled_) {
    printEnergies();
  }
  if (isCheckAbortCriteraSlipRateEnabled_) {
    checkAbortCriterion(minTimeSinceSlipRateBelowThreshold_, "All slip rates are");
  }
  if (isCheckAbortCriteraMomentRateEnabled_) {
    checkAbortCriterion(minTimeSinceMomentRateBelowThreshold_, "The seismic moment rate is");
  }

  if (isFileOutputEnabled_) {
    writeEnergies(time);
  }
  ++outputId_;
  logInfo() << "Writing energy output at time" << time << "Done.";
}

void EnergyOutput::simulationStart(std::optional<double> checkpointTime) {
  if (isFileOutputEnabled_) {
    out_.open(outputFileName_);
    out_ << std::scientific;
    out_ << std::setprecision(std::numeric_limits<double>::max_digits10);
    writeHeader();
  }
  syncPoint(checkpointTime.value_or(0));
}

EnergyOutput::~EnergyOutput() = default;

void EnergyOutput::computeDynamicRuptureEnergies() {
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    double& totalFrictionalWork = energiesStorage_.totalFrictionalWork(sim);
    double& staticFrictionalWork = energiesStorage_.staticFrictionalWork(sim);
    double& seismicMoment = energiesStorage_.seismicMoment(sim);
    double& potency = energiesStorage_.potency(sim);
    minTimeSinceSlipRateBelowThreshold_[sim] = std::numeric_limits<double>::infinity();

    for (const auto& layer : drStorage_->leaves()) {

      real* const* timeDofsPlus = layer.var<DynamicRupture::TimeDerivativePlus>();
      real* const* timeDofsMinus = layer.var<DynamicRupture::TimeDerivativeMinus>();

      const auto* godunovData = layer.var<DynamicRupture::GodunovData>();
      const auto* faceInformation = layer.var<DynamicRupture::FaceInformation>();
      const auto* drEnergyOutput = layer.var<DynamicRupture::DREnergyOutputVar>();
      const auto* waveSpeedsPlus = layer.var<DynamicRupture::WaveSpeedsPlus>();
      const auto* waveSpeedsMinus = layer.var<DynamicRupture::WaveSpeedsMinus>();
      const auto layerSize = layer.size();

#if !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(                                                                \
        + : totalFrictionalWork, staticFrictionalWork, seismicMoment, potency) default(none)       \
    shared(layerSize,                                                                              \
               drEnergyOutput,                                                                     \
               faceInformation,                                                                    \
               timeDofsMinus,                                                                      \
               timeDofsPlus,                                                                       \
               godunovData,                                                                        \
               waveSpeedsPlus,                                                                     \
               waveSpeedsMinus,                                                                    \
               sim)
#endif
      for (std::size_t i = 0; i < layerSize; ++i) {
        if (faceInformation[i].plusSideOnThisRank) {

          for (std::size_t j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints; ++j) {
            totalFrictionalWork +=
                drEnergyOutput[i].frictionalEnergy[j * seissol::multisim::NumSimulations + sim];
          }
          staticFrictionalWork += computeStaticWork(timeDofsPlus[i],
                                                    timeDofsMinus[i],
                                                    faceInformation[i],
                                                    godunovData[i],
                                                    drEnergyOutput[i].slip,
                                                    global_)[sim];

          const double muPlus = waveSpeedsPlus[i].density * waveSpeedsPlus[i].sWaveVelocity *
                                waveSpeedsPlus[i].sWaveVelocity;
          const double muMinus = waveSpeedsMinus[i].density * waveSpeedsMinus[i].sWaveVelocity *
                                 waveSpeedsMinus[i].sWaveVelocity;
          const double mu = 2.0 * muPlus * muMinus / (muPlus + muMinus);
          double potencyIncrease = 0.0;
          for (std::size_t k = 0; k < seissol::dr::misc::NumBoundaryGaussPoints; ++k) {
            potencyIncrease +=
                drEnergyOutput[i].accumulatedSlip[k * seissol::multisim::NumSimulations + sim];
          }
          potencyIncrease *=
              0.5 * godunovData[i].doubledSurfaceArea / seissol::dr::misc::NumBoundaryGaussPoints;
          potency += potencyIncrease;
          seismicMoment += potencyIncrease * mu;
        }
      }
      double localMin = std::numeric_limits<double>::infinity();
#if !NVHPC_AVOID_OMP
#pragma omp parallel for reduction(min : localMin) default(none)                                   \
    shared(layerSize, drEnergyOutput, faceInformation, sim)
#endif
      for (std::size_t i = 0; i < layerSize; ++i) {
        if (faceInformation[i].plusSideOnThisRank) {
          for (std::size_t j = 0; j < seissol::dr::misc::NumBoundaryGaussPoints; ++j) {
            localMin = std::min(
                static_cast<double>(
                    drEnergyOutput[i].timeSinceSlipRateBelowThreshold
                        [static_cast<size_t>(j * seissol::multisim::NumSimulations) + sim]),
                localMin);
          }
        }
      }
      minTimeSinceSlipRateBelowThreshold_[sim] =
          std::min(localMin, minTimeSinceSlipRateBelowThreshold_[sim]);
    }
  }
}

void EnergyOutput::computeVolumeEnergies() {
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    auto& totalGravitationalEnergyLocal = energiesStorage_.gravitationalEnergy(sim);
    auto& totalAcousticEnergyLocal = energiesStorage_.acousticEnergy(sim);
    auto& totalAcousticKineticEnergyLocal = energiesStorage_.acousticKineticEnergy(sim);
    auto& totalElasticEnergyLocal = energiesStorage_.elasticEnergy(sim);
    auto& totalElasticKineticEnergyLocal = energiesStorage_.elasticKineticEnergy(sim);
    auto& totalPlasticMoment = energiesStorage_.plasticMoment(sim);
    auto& totalMomentumXLocal = energiesStorage_.totalMomentumX(sim);
    auto& totalMomentumYLocal = energiesStorage_.totalMomentumY(sim);
    auto& totalMomentumZLocal = energiesStorage_.totalMomentumZ(sim);

    const std::unordered_map<std::string, double&> energyRefMap{
        {"acoustic-energy", totalAcousticEnergyLocal},
        {"acoustic-kinetic-energy", totalAcousticKineticEnergyLocal},
        {"elastic-energy", totalElasticEnergyLocal},
        {"elastic-kinetic-energy", totalElasticKineticEnergyLocal},
        {"momentum-x", totalMomentumXLocal},
        {"momentum-y", totalMomentumYLocal},
        {"momentum-z", totalMomentumZLocal}};

    const std::vector<Element>& elements = meshReader_->getElements();
    const std::vector<Vertex>& vertices = meshReader_->getVertices();

    const auto g = seissolInstance_.getGravitationSetup().acceleration;

    constexpr auto QuadPolyDegree = ConvergenceOrder + 1;
    constexpr auto NumQuadraturePointsTet = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

    double quadraturePointsTet[NumQuadraturePointsTet][3]{};
    double quadratureWeightsTet[NumQuadraturePointsTet]{};
    seissol::quadrature::TetrahedronQuadrature(
        quadraturePointsTet, quadratureWeightsTet, QuadPolyDegree);

    constexpr auto NumQuadraturePointsTri = QuadPolyDegree * QuadPolyDegree;
    double quadraturePointsTri[NumQuadraturePointsTri][2]{};
    double quadratureWeightsTri[NumQuadraturePointsTri]{};
    seissol::quadrature::TriangleQuadrature(
        quadraturePointsTri, quadratureWeightsTri, QuadPolyDegree);

    // Note: Default(none) is not possible, clang requires data sharing attribute for g, gcc forbids
    // it
    for (const auto& layer : ltsStorage_->leaves(Ghost)) {
      const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
      const auto* cellInformationData = layer.var<LTS::CellInformation>();
      const auto* faceDisplacementsData = layer.var<LTS::FaceDisplacements>();
      const auto* materialData = layer.var<LTS::MaterialData>();
      const auto* boundaryMappingData = layer.var<LTS::BoundaryMapping>();
      const auto* pstrainData = layer.var<LTS::PStrain>();
      const auto* dofsData = layer.var<LTS::Dofs>();
      const auto* energyData = layer.var<LTS::EnergyData>();

      constexpr auto EnergyCount = model::EnergyCompute<model::MaterialT>::EnergyCount;
      double energyValues[EnergyCount]{};

#if !NVHPC_AVOID_OMP
#pragma omp parallel for schedule(static)                                                          \
    reduction(+ : totalGravitationalEnergyLocal, energyValues[ : EnergyCount], totalPlasticMoment) \
    shared(elements, vertices, global_)
#endif
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        if (secondaryInformation[cell].duplicate > 0) {
          // skip duplicate cells
          continue;
        }
        const auto elementId = secondaryInformation[cell].meshId;
        const double volume = MeshTools::volume(elements[elementId], vertices);

        // NOLINTNEXTLINE
        const auto& material = materialData[cell];
        const auto& cellInformation = cellInformationData[cell];
        const auto& faceDisplacements = faceDisplacementsData[cell];

        // Needed to weight the integral.
        const auto jacobiDet = 6 * volume;

        alignas(Alignment) real linData[tensor::massLPR::size()];
        auto lin = init::massLPR::view::create(linData);
        // Evaluate numerical solution at quad. nodes
        kernel::massLP krnl;
        krnl.M3 = init::M3::Values;
        krnl.massLPR = linData;
        krnl.Q = dofsData[cell];
        krnl.execute();

        alignas(Alignment) real quadData[tensor::massSPR::size()];
        auto quad = init::massSPR::view::create(quadData);
        // Evaluate numerical solution at quad. nodes
        kernel::massSP krnl2;
        krnl2.M3 = init::M3::Values;
        krnl2.massSPR = quadData;
        krnl2.Q = dofsData[cell];
        krnl2.execute();

        auto linSub = multisim::simtensor(lin, sim);
        auto quadSub = multisim::simtensor(quad, sim);

        // assume _constant_ material over a cell (will need adjustments for e.g. #1297)

        const auto localValues = model::EnergyCompute<model::MaterialT>::computeEnergies(
            material, energyData[cell], linSub, quadSub);

        for (std::size_t i = 0; i < localValues.size(); ++i) {
          energyValues[i] += jacobiDet * localValues[i];
        }

        constexpr auto UIdx = model::MaterialT::TractionQuantities;

        const auto& boundaryMappings = boundaryMappingData[cell];
        // Compute gravitational energy
        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          if (cellInformation.faceTypes[face] != FaceType::FreeSurfaceGravity) {
            continue;
          }

          // Displacements are stored in face-aligned coordinate system.
          // We need to rotate it to the global coordinate system.
          const auto& boundaryMapping = boundaryMappings[face];
          auto tinv = init::Tinv::view::create(boundaryMapping.dataTinv);
          alignas(Alignment)
              real rotateDisplacementToFaceNormalData[init::displacementRotationMatrix::Size];

          auto rotateDisplacementToFaceNormal =
              init::displacementRotationMatrix::view::create(rotateDisplacementToFaceNormalData);
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              rotateDisplacementToFaceNormal(i, j) = tinv(i + UIdx, j + UIdx);
            }
          }

          const auto* curFaceDisplacementsData = faceDisplacements[face];

          // See for example (Saito, Tsunami generation and propagation, 2019) section 3.2.3 for
          // derivation.

          alignas(Alignment) std::array<real, tensor::quadFaceDisplacements::Size>
              quadFaceDisplacements{};
          {
            seissol::kernel::quadFaceDisplacementsCompute evalKrnl;
            evalKrnl.rotatedFaceDisplacement = curFaceDisplacementsData;
            evalKrnl.M2 = init::M2::Values;
            evalKrnl.MV2nTo2m = nodal::init::MV2nTo2m::Values;
            evalKrnl.quadFaceDisplacements = quadFaceDisplacements.data();
            evalKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
            evalKrnl.execute();
          }

          const auto quadFaceDisplacementsViewFused =
              init::quadFaceDisplacements::view::create(quadFaceDisplacements.data());
          const auto quadFaceDisplacementsView =
              multisim::simtensor(quadFaceDisplacementsViewFused, sim);

          const auto surface = MeshTools::surface(elements[elementId], face, vertices);
          const auto rho = material.getDensity();

          // contains an elided 0.5 * 2.0 (1/2 due to energy; 2 due to surface)
          totalGravitationalEnergyLocal += rho * g * surface * quadFaceDisplacementsView(0);
        }

        if (isPlasticityEnabled_) {
          // plastic moment
          const real* pstrainCell = pstrainData[cell];
          const double mu = material.getMuBar();

          // integrating over all collocation points suffices
          const real* __restrict qEta = &pstrainCell[tensor::QStressNodal::size()];

          alignas(Alignment) real qEtaQuad[tensor::QEtaNodalProject::size()]{};

          kernel::plProject krnl;
          set_evalAtQP(krnl, global_->evalAtQPMatrix);
          set_vInv(krnl, global_->vandermondeMatrixInverse);
          krnl.QEtaNodal = qEta;
          krnl.QEtaNodalProject = qEtaQuad;
          krnl.execute();

          double pMoment = 0;

#pragma omp simd reduction(+ : pMoment)
          for (size_t qp = 0; qp < NumQuadraturePointsTet; ++qp) {
            pMoment += quadratureWeightsTet[qp] * qEtaQuad[qp];
          }

          totalPlasticMoment += mu * jacobiDet * pMoment;
        }
      }

      for (std::size_t i = 0; i < model::EnergyCompute<model::MaterialT>::EnergyCount; ++i) {
        energyRefMap.at(model::EnergyCompute<model::MaterialT>::Energies[i]) += energyValues[i];
      }
    }
  }
}

void EnergyOutput::computeEnergies() {
  energiesStorage_.energies.fill(0.0);
  if (shouldComputeVolumeEnergies()) {
    computeVolumeEnergies();
  }
  computeDynamicRuptureEnergies();
}

void EnergyOutput::reduceEnergies() {
  const auto& comm = Mpi::mpi.comm();
  MPI_Allreduce(MPI_IN_PLACE,
                energiesStorage_.energies.data(),
                static_cast<int>(energiesStorage_.energies.size()),
                MPI_DOUBLE,
                MPI_SUM,
                comm);
}

void EnergyOutput::reduceMinTimeSinceSlipRateBelowThreshold() {
  const auto& comm = Mpi::mpi.comm();
  MPI_Allreduce(MPI_IN_PLACE,
                minTimeSinceSlipRateBelowThreshold_.data(),
                static_cast<int>(minTimeSinceSlipRateBelowThreshold_.size()),
                Mpi::castToMpiType<double>(),
                MPI_MIN,
                comm);
}

void EnergyOutput::printEnergies() {
  const auto outputPrecision =
      seissolInstance_.getSeisSolParameters().output.energyParameters.terminalPrecision;

  const auto shouldPrint = [](double thresholdValue) { return std::abs(thresholdValue) > 1.e-20; };
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    const std::string fusedPrefix =
        multisim::MultisimEnabled ? "[" + std::to_string(sim) + "]" : "";
    const auto totalAcousticEnergy =
        energiesStorage_.acousticKineticEnergy(sim) + energiesStorage_.acousticEnergy(sim);
    const auto totalElasticEnergy =
        energiesStorage_.elasticKineticEnergy(sim) + energiesStorage_.elasticEnergy(sim);
    const auto ratioElasticKinematic =
        100.0 * energiesStorage_.elasticKineticEnergy(sim) / totalElasticEnergy;
    const auto ratioElasticPotential =
        100.0 * energiesStorage_.elasticEnergy(sim) / totalElasticEnergy;
    const auto ratioAcousticKinematic =
        100.0 * energiesStorage_.acousticKineticEnergy(sim) / totalAcousticEnergy;
    const auto ratioAcousticPotential =
        100.0 * energiesStorage_.acousticEnergy(sim) / totalAcousticEnergy;
    const auto totalFrictionalWork = energiesStorage_.totalFrictionalWork(sim);
    const auto staticFrictionalWork = energiesStorage_.staticFrictionalWork(sim);
    const auto radiatedEnergy = totalFrictionalWork - staticFrictionalWork;
    const auto ratioFrictionalStatic = 100.0 * staticFrictionalWork / totalFrictionalWork;
    const auto ratioFrictionalRadiated = 100.0 * radiatedEnergy / totalFrictionalWork;
    const auto ratioPlasticMoment =
        100.0 * energiesStorage_.plasticMoment(sim) /
        (energiesStorage_.plasticMoment(sim) + energiesStorage_.seismicMoment(sim));
    const auto totalMomentumX = energiesStorage_.totalMomentumX(sim);
    const auto totalMomentumY = energiesStorage_.totalMomentumY(sim);
    const auto totalMomentumZ = energiesStorage_.totalMomentumZ(sim);
    if (shouldComputeVolumeEnergies()) {
      if (shouldPrint(totalElasticEnergy)) {
        logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                  << "Elastic energy (total, % kinematic, % potential): "
                  << UnitEnergy.formatScientific(totalElasticEnergy, {}, outputPrecision).c_str()
                  << "," << ratioElasticKinematic << "% ," << ratioElasticPotential << "%";
      }
      if (shouldPrint(totalAcousticEnergy)) {
        logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                  << "Acoustic energy (total, % kinematic, % potential): "
                  << UnitEnergy.formatScientific(totalAcousticEnergy, {}, outputPrecision).c_str()
                  << "," << ratioAcousticKinematic << "% ," << ratioAcousticPotential << "%";
      }
      if (shouldPrint(energiesStorage_.gravitationalEnergy(sim))) {
        logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                  << "Gravitational energy:"
                  << UnitEnergy
                         .formatScientific(
                             energiesStorage_.gravitationalEnergy(sim), {}, outputPrecision)
                         .c_str();
      }
      if (shouldPrint(energiesStorage_.plasticMoment(sim))) {
        logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                  << "Plastic moment (value, equivalent Mw, % total moment):"
                  << UnitMoment
                         .formatScientific(energiesStorage_.plasticMoment(sim), {}, outputPrecision)
                         .c_str()
                  << "," << 2.0 / 3.0 * std::log10(energiesStorage_.plasticMoment(sim)) - 6.07
                  << "," << ratioPlasticMoment << "%";
      }
    } else {
      logInfo() << "Volume energies skipped at this step";
    }
    logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
              << " Total momentum (X, Y, Z):"
              << UnitMomentum.formatScientific(totalMomentumX, {}, outputPrecision).c_str() << ","
              << UnitMomentum.formatScientific(totalMomentumY, {}, outputPrecision).c_str() << ","
              << UnitMomentum.formatScientific(totalMomentumZ, {}, outputPrecision).c_str();
    if (shouldPrint(totalFrictionalWork)) {
      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << "Frictional work (total, % static, % radiated): "
                << UnitEnergy.formatScientific(totalFrictionalWork, {}, outputPrecision).c_str()
                << "," << ratioFrictionalStatic << "% ," << ratioFrictionalRadiated << "%";
      logInfo() << std::setprecision(outputPrecision) << fusedPrefix.c_str()
                << "Seismic moment (without plasticity):"
                << UnitMoment
                       .formatScientific(energiesStorage_.seismicMoment(sim), {}, outputPrecision)
                       .c_str()
                << ", Mw:" << 2.0 / 3.0 * std::log10(energiesStorage_.seismicMoment(sim)) - 6.07;
    }
    if (!std::isfinite(totalElasticEnergy + totalAcousticEnergy)) {
      logError() << fusedPrefix << " Detected Inf/NaN in energies. Aborting.";
    }
  }
}

void EnergyOutput::checkAbortCriterion(
    const std::array<double, multisim::NumSimulations>& timeSinceThreshold,
    const std::string& prefixMessage) {
  size_t abortCount = 0;
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    if ((timeSinceThreshold[sim] > 0) and
        (timeSinceThreshold[sim] < std::numeric_limits<double>::infinity())) {
      if (static_cast<double>(timeSinceThreshold[sim]) < terminatorMaxTimePostRupture_) {
        logInfo() << prefixMessage.c_str() << "below threshold since" << timeSinceThreshold[sim]
                  << "s; in simulation: " << sim
                  << "(lower than the abort criteria: " << terminatorMaxTimePostRupture_ << "s)";
      } else {
        logInfo() << prefixMessage.c_str() << "below threshold since" << timeSinceThreshold[sim]
                  << "s; in simulation: " << sim
                  << "(greater than the abort criteria: " << terminatorMaxTimePostRupture_ << "s)";
        ++abortCount;
      }
    }
  }

  bool abort = abortCount == multisim::NumSimulations;
  const auto& comm = Mpi::mpi.comm();
  MPI_Bcast(reinterpret_cast<void*>(&abort), 1, MPI_CXX_BOOL, 0, comm);
  if (abort) {
    seissolInstance_.simulator().abort();
  }
}

void EnergyOutput::writeHeader() {
  out_ << "time,variable,simulation_index,measurement" << std::endl;
}

void EnergyOutput::writeEnergies(double time) {
  for (size_t sim = 0; sim < multisim::NumSimulations; sim++) {
    const std::string fusedSuffix = std::to_string(sim);
    if (shouldComputeVolumeEnergies()) {
      out_ << time << ",gravitational_energy," << fusedSuffix << ","
           << energiesStorage_.gravitationalEnergy(sim) << "\n"
           << time << ",acoustic_energy," << fusedSuffix << ","
           << energiesStorage_.acousticEnergy(sim) << "\n"
           << time << ",acoustic_kinetic_energy," << fusedSuffix << ","
           << energiesStorage_.acousticKineticEnergy(sim) << "\n"
           << time << ",elastic_energy," << fusedSuffix << ","
           << energiesStorage_.elasticEnergy(sim) << "\n"
           << time << ",elastic_kinetic_energy," << fusedSuffix << ","
           << energiesStorage_.elasticKineticEnergy(sim) << "\n"
           << time << ",plastic_moment," << fusedSuffix << ","
           << energiesStorage_.plasticMoment(sim) << "\n"
           << time << ",momentumX," << fusedSuffix << "," << energiesStorage_.totalMomentumX(sim)
           << "\n"
           << time << ",momentumY," << fusedSuffix << "," << energiesStorage_.totalMomentumY(sim)
           << "\n"
           << time << ",momentumZ," << fusedSuffix << "," << energiesStorage_.totalMomentumZ(sim)
           << "\n";
    }
    out_ << time << ",total_frictional_work," << fusedSuffix << ","
         << energiesStorage_.totalFrictionalWork(sim) << "\n"
         << time << ",static_frictional_work," << fusedSuffix << ","
         << energiesStorage_.staticFrictionalWork(sim) << "\n"
         << time << ",seismic_moment," << fusedSuffix << "," << energiesStorage_.seismicMoment(sim)
         << "\n"
         << time << ",potency," << fusedSuffix << "," << energiesStorage_.potency(sim) << "\n"
         << time << ",plastic_moment," << fusedSuffix << "," << energiesStorage_.plasticMoment(sim)
         << std::endl;
  }
}

bool EnergyOutput::shouldComputeVolumeEnergies() const {
  return outputId_ % computeVolumeEnergiesEveryOutput_ == 0;
}

} // namespace seissol::writer
