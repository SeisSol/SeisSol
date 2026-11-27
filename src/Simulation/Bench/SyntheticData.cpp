// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SyntheticData.h"

#include "Alignment.h"
#include "Config.h"
#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Typedefs.h"
#include "Core/Memory/CellLocalInformation.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/GlobalData.h"
#include "Core/Memory/LtsSetup.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/SimulationSettings.h"
#include "Core/Memory/Tree/Colormap.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Parallel/OpenMP.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Physics/Impedances.h"
#include "Features/DynamicRupture/Setup/Factory.h"
#include "Features/DynamicRupture/Setup/Initializer/BaseDRInitializer.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "Features/DynamicRupture/Typedefs.h"
#include "Features/Plasticity/Solver/Storage.h"
#include "GeneratedCode/tensor.h"
#include "Simulation/Bench/LoopBench.h"
#include "Simulation/Bench/Proxy.h"
#include "Simulation/InitProcedure/Internal/Recording.h"
#include "WavePropagation/Physics/CommonDatastructures.h"
#include "WavePropagation/Solver/SolverSelector.h"
#include "WavePropagation/Solver/Storage.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#ifdef ACL_DEVICE
#include "Simulation/InitProcedure/Internal/Scratchpads.h"

#include <Device/device.h>
#endif

namespace seissol::proxy {

namespace {

/// Fills `buffer` with arbitrary finite values that are not zero.
template <typename RealT>
void fillWithStuff(RealT* buffer, std::size_t count) {
  // No real point for these numbers. Should be just something != 0 and != NaN and != Inf
  const auto stuff = [](std::size_t n) {
    return static_cast<RealT>((214013.0 * static_cast<double>(n) + 2531011.0) / 16777216.0);
  };
#pragma omp parallel for schedule(static)
  for (std::size_t n = 0; n < count; ++n) {
    buffer[n] = stuff(n);
  }
}

/// the material on both sides of the fault faces: that of the SCEC benchmarks TPV5, TPV101, ...
constexpr model::IsotropicWaveSpeeds FaultMaterial{
    .density = 2670.0, .pWaveVelocity = 6000.0, .sWaveVelocity = 3464.0};

/// the normal stress on the fault faces, a compression
constexpr double NormalStress = -120.0e6;

/// The parameters of the parameter file for the friction law `frictionLaw`, from a benchmark of it:
/// TPV101 for the aging law and its siblings, TPV104 for the fast velocity weakening, TPV105 for
/// thermal pressurization.
auto drParametersOf(FrictionLaw frictionLaw)
    -> std::shared_ptr<initializer::parameters::DRParameters> {
  auto parameters = std::make_shared<initializer::parameters::DRParameters>();
  parameters->frictionLawType = static_cast<dr::misc::FrictionLawType>(frictionLaw.type);
  parameters->isThermalPressureOn = frictionLaw.thermalPressurization;
  switch (frictionLaw.type) {
  case 3:
  case 4:
  case 7:
    parameters->rsF0 = 0.6;
    parameters->rsSr0 = 1e-6;
    parameters->rsB = 0.012;
    parameters->rsInitialSlipRate1 = 1e-16;
    break;
  case 103:
    parameters->rsF0 = 0.6;
    parameters->rsSr0 = 1e-6;
    parameters->rsB = 0.014;
    parameters->muW = 0.2;
    parameters->rsInitialSlipRate1 = 1e-16;
    break;
  case 6:
    parameters->vStar = 0.1;
    parameters->prakashLength = 0.01;
    break;
  default:
    break;
  }
  if (frictionLaw.thermalPressurization) {
    parameters->thermalDiffusivity = 1e-6;
    parameters->heatCapacity = 2.7e6;
    parameters->undrainedTPResponse = 0.1e6;
    parameters->initialTemperature = 483.15;
    parameters->initialPressure = 0;
  }
  return parameters;
}

/// The shear traction on the fault faces under the friction law `type`: from a benchmark of it, as
/// for drParametersOf(), above the strength for the linear slip weakening (TPV5 in its nucleation
/// patch); none for the imposed slip rates.
auto shearStressOf(std::uint32_t type) -> double {
  switch (type) {
  case 6:
  case 16:
  case 1058:
    return 81.6e6;
  case 3:
  case 4:
    return 75.0e6;
  case 7:
    return 70.8e6;
  case 103:
    return 40.0e6;
  case 33:
  case 34:
  case 35:
    return 0;
  default:
    return 70.0e6;
  }
}

/// The rise of the shear stress after the initial state of the laws of rate and state: their
/// nucleation in the benchmark (TPV101, TPV104), which overstresses all faces here, not a patch.
auto nucleationStressOf(std::uint32_t type) -> double {
  switch (type) {
  case 3:
  case 4:
  case 7:
    return 25.0e6;
  case 103:
    return 45.0e6;
  default:
    return 0;
  }
}

/// The parameters of the fault files for the friction law `type`, by their names there (see
/// BaseDRInitializer::lawParameters()): from a benchmark of it, as for drParametersOf().
auto faultParametersOf(std::uint32_t type) -> std::unordered_map<std::string, double> {
  std::unordered_map<std::string, double> values{
      // linear slip weakening (TPV5)
      {"mu_s", 0.677},
      {"mu_d", 0.525},
      {"d_c", 0.4},
      {"cohesion", 0.0},
      // thermal pressurization (TPV105)
      {"tp_halfWidthShearZone", 0.02},
      {"tp_hydraulicDiffusivity", 4e-4},
      // the source time functions of the imposed slip rates
      {"tau_S", 0.05},
      {"tau_R", 1.0},
      {"rupture_rise_time", 1.0},
  };
  // rate and state
  if (type == 103) {
    values.insert({{"rs_a", 0.01}, {"rs_sl0", 0.4}, {"rs_srW", 0.1}});
  } else if (type == 7) {
    values.insert({{"rs_a", 0.008}, {"rs_sl0", 1e-7}});
  } else {
    values.insert({{"rs_a", 0.008}, {"rs_sl0", 0.02}});
  }
  return values;
}

} // namespace

template <typename Cfg>
SyntheticData<Cfg>::SyntheticData(const BenchSetup& setup)
    : cellCount_(setup.cells), dynamicRupture_(setup.dynamicRupture), plasticity_(setup.plasticity),
      layerId_(initializer::LayerIdentifier(HaloType::Interior, configIdOf<Cfg>(), 0)) {
  initGlobalData();
  initCells(dynamicRupture_ ? FaceType::DynamicRupture : FaceType::Regular);
  if (plasticity_ != CellPlasticity::None) {
    initPlasticity();
  }
  if (dynamicRupture_) {
    initFaults(setup.frictionLaw);
  }
  initDevice();
}

template <typename Cfg>
void SyntheticData<Cfg>::initGlobalData() {
  initializer::GlobalDataInitializerOnHost::init<Cfg>(
      globalData_.template onHost<Cfg>(), allocator_, memory::Memkind::Standard);
  if constexpr (isDeviceOn()) {
    initializer::GlobalDataInitializerOnDevice::init<Cfg>(
        globalData_.template onDevice<Cfg>(), allocator_, memory::Memkind::DeviceGlobalMemory);
  }
}

template <typename Cfg>
void SyntheticData<Cfg>::initCells(FaceType faceType) {
  const initializer::LTSColorMap map(initializer::EnumLayer<HaloType>({HaloType::Interior}),
                                     initializer::EnumLayer<std::size_t>({0}),
                                     initializer::EnumLayer<ConfigId>({configIdOf<Cfg>()}));

  const bool plasticity = plasticity_ != CellPlasticity::None;
  const SimulationSettings settings(plasticity, false, Cfg::NumSimulations, Cfg::NumSimulations);
  LTS::addTo(ltsStorage_);
  LTS::addWavePropagation(ltsStorage_, settings);
  LTSPlasticity::addTo(ltsStorage_, plasticity);
  ltsStorage_.setLayerCount(map);
  ltsStorage_.fixate();

  LTS::Layer& layer = cells();
  layer.setNumberOfCells(cellCount_);
  layer.setEntrySize<LTS::Buffers>(sizeof(real) * kernels::SolverOf<Cfg>::IntegralsSize *
                                   cellCount_);

  ltsStorage_.allocateVariables();
  ltsStorage_.touchVariables();
  ltsStorage_.allocateBuckets();

  real(*dofs)[tensor::Q<Cfg>::size()] = layer.var<LTS::Dofs>(Cfg());
  real** buffers = layer.var<LTS::StepIntegrals>(Cfg());
  real** derivatives = layer.var<LTS::Derivatives>(Cfg());
  auto* faceNeighbors = layer.var<LTS::FaceNeighbors>();
  auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
  real* bucket =
      static_cast<real*>(layer.var<LTS::Buffers>(Cfg(), initializer::AllocationPlace::Host));

  // every cell keeps the time integral of its step in its buffer, which its neighbors read
  LtsSetup gts;
  gts.setHasBuffer(true, BufferType::StepIntegrals);
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    gts.setNeighborBuffer(face, BufferType::StepIntegrals).setNeighborGTSRelation(face, true);
  }

  std::mt19937 rng(layer.size());
  std::uniform_int_distribution<unsigned> sideDist(0, 3);
  std::uniform_int_distribution<std::size_t> cellDist(0, layer.size() - 1);

  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    buffers[cell] = bucket + cell * kernels::SolverOf<Cfg>::IntegralsSize;
    derivatives[cell] = nullptr;

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      cellInformation[cell].faceTypes[face] = faceType;
      cellInformation[cell].faceRelations[face][0] = sideDist(rng);
      cellInformation[cell].faceRelations[face][1] = 0;
      cellInformation[cell].neighborConfigIds[face] = configIdOf<Cfg>();

      const auto neighbor = cellDist(rng);
      secondaryInformation[cell].faceNeighbors[face].global = neighbor;
      secondaryInformation[cell].faceNeighbors[face].color = 0;
      secondaryInformation[cell].faceNeighbors[face].cell = neighbor;
    }
    cellInformation[cell].ltsSetup = gts;
  }

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (faceType == FaceType::Regular) {
        const auto neighbor = secondaryInformation[cell].faceNeighbors[face].cell;
        faceNeighbors[cell][face] = buffers[neighbor];
      } else {
        faceNeighbors[cell][face] = nullptr;
      }
    }
  }

  // the same addresses on the device (without one, these variables stay empty)
  if constexpr (isDeviceOn()) {
    real** buffersDevice = layer.var<LTS::StepIntegralsDevice>(Cfg());
    real** derivativesDevice = layer.var<LTS::DerivativesDevice>(Cfg());
    auto* faceNeighborsDevice = layer.var<LTS::FaceNeighborsDevice>();
    real* bucketDevice =
        static_cast<real*>(layer.var<LTS::Buffers>(Cfg(), initializer::AllocationPlace::Device));
#pragma omp parallel for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      buffersDevice[cell] = bucketDevice + cell * kernels::SolverOf<Cfg>::IntegralsSize;
      derivativesDevice[cell] = nullptr;
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (faceType == FaceType::Regular) {
          const auto neighbor = secondaryInformation[cell].faceNeighbors[face].cell;
          faceNeighborsDevice[cell][face] =
              bucketDevice + neighbor * kernels::SolverOf<Cfg>::IntegralsSize;
        }
      }
    }
  }

  fillWithStuff(reinterpret_cast<real*>(dofs), tensor::Q<Cfg>::size() * layer.size());
  fillWithStuff(bucket, kernels::SolverOf<Cfg>::IntegralsSize * layer.size());
  // The local and the neighboring integration data (the star matrices, the flux solvers, ...)
  // stay zero: the kernels compute with them as with any other values, and the cells keep the
  // values of their unknowns from step to step. With arbitrary values there, these grew by about
  // 24 orders of magnitude in every step, and were not finite after a few (in single precision,
  // after the first); the plasticity then found nothing to yield.

  if constexpr (Cfg::Solver == SolverType::STP) {
    auto* localIntegration = layer.var<LTS::LocalIntegration>(Cfg());
#pragma omp parallel for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      localIntegration[cell].specific.typicalTimeStepWidth = Timestep;
    }
  }
}

template <typename Cfg>
void SyntheticData<Cfg>::initPlasticity() {
  LTS::Layer& layer = cells();
  auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* plasticity = layer.var<LTSPlasticity::Plasticity>(Cfg());

  // A cell without cohesion and friction yields wherever its stress has a deviator, which the
  // arbitrary values of the cells give it; with the slow relaxation (see RelaxationTime), it
  // stays beyond the yield surface. Otherwise, the cohesion lies far above any stress of the
  // cells.
  const auto cohesion =
      plasticity_ == CellPlasticity::Yield ? static_cast<real>(0) : static_cast<real>(1e20);
  // 1 / (2 mu), with the shear modulus of crustal rock (32 GPa)
  const auto muFactor = static_cast<real>(1.0 / (2.0 * 3.2e10));

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    cellInformation[cell].plasticityEnabled = true;
    auto& data = plasticity[cell];
    std::fill(std::begin(data.initialLoading), std::end(data.initialLoading), 0);
    std::fill(std::begin(data.cohesionTimesCosAngularFriction),
              std::end(data.cohesionTimesCosAngularFriction),
              cohesion);
    std::fill(std::begin(data.sinAngularFriction), std::end(data.sinAngularFriction), 0);
    data.mufactor = muFactor;
  }
}

template <typename Cfg>
void SyntheticData<Cfg>::initFaults(FrictionLaw frictionLaw) {
  // the variables and the solvers of the friction law, as a run makes them, without the mesh
  drParameters_ = drParametersOf(frictionLaw);
  const auto factory = dr::factory::getFactory(drParameters_);
  auto solver = factory->solver();
  faultVariables_ = std::move(solver.storage);
  frictionLaw_ = std::move(solver.frictionLaw);
  frictionLawDevice_ = std::move(solver.frictionLawDevice);

  faultVariables_->addTo(drStorage_);
  drStorage_.setLayerCount(ltsStorage_.getColorMap());
  drStorage_.fixate();

  DynamicRupture::Layer& faults = faces();
  faults.setNumberOfCells(Cell::NumFaces * cellCount_);

  drStorage_.allocateVariables();
  drStorage_.touchVariables();

  constexpr auto DerivativesSize = kernels::SolverOf<Cfg>::DerivativesSize;
  const auto derivativesBytes = cellCount_ * DerivativesSize * sizeof(real);
  faceDerivatives_ = static_cast<real*>(
      allocator_.allocateMemory(derivativesBytes, PagesizeHeap, memory::Memkind::Standard));

#pragma omp parallel
  {
    std::mt19937 rng(cellCount_ + OpenMP::threadId());
    std::uniform_real_distribution<real> urd;
#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < cellCount_; ++cell) {
      for (std::size_t i = 0; i < DerivativesSize; ++i) {
        faceDerivatives_[cell * DerivativesSize + i] = urd(rng);
      }
    }
  }

#ifdef ACL_DEVICE
  faceDerivativesOnDevice_ = static_cast<real*>(allocator_.allocateMemory(
      derivativesBytes, PagesizeHeap, memory::Memkind::DeviceGlobalMemory));
  device::DeviceInstance::instance().api().copyTo(
      faceDerivativesOnDevice_, faceDerivatives_, derivativesBytes);
#endif

  LTS::Layer& layer = cells();
  auto* drMapping = layer.var<LTS::DRMapping>(Cfg());

  real(*imposedStatePlus)[tensor::QInterpolated<Cfg>::size()] =
      faults.var<DynamicRupture::ImposedStatePlus>(Cfg());
  real(*fluxSolverPlus)[tensor::fluxSolver<Cfg>::size()] =
      faults.var<DynamicRupture::FluxSolverPlus>(Cfg());
  real** timeDerivativePlus = faults.var<DynamicRupture::TimeDerivativePlus>(Cfg());
  real** timeDerivativeMinus = faults.var<DynamicRupture::TimeDerivativeMinus>(Cfg());
  auto* faceInformation = faults.var<DynamicRupture::FaceInformation>();

  std::mt19937 rng(cellCount_);
  std::uniform_int_distribution<unsigned> sideDist(0, 3);
  std::uniform_int_distribution<unsigned> orientationDist(0, 1);
  std::uniform_int_distribution<std::size_t> faceDist(0, faults.size() - 1);
  std::uniform_int_distribution<std::size_t> cellDist(0, cellCount_ - 1);

  // each face of a cell takes the flux of a random fault face; as in the solver, the mapping of a
  // face addresses that face
  for (std::size_t cell = 0; cell < cellCount_; ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto orientation = orientationDist(rng);
      const auto faultFace = faceDist(rng);
      auto& mapping = drMapping[cell][face];
      mapping.side = static_cast<std::int8_t>(face);
      mapping.faceRelation = static_cast<std::int8_t>(orientation);
      mapping.godunov = imposedStatePlus[faultFace];
      mapping.fluxSolver = fluxSolverPlus[faultFace];
    }
  }

  for (std::size_t face = 0; face < faults.size(); ++face) {
    const auto plusCell = cellDist(rng);
    const auto minusCell = cellDist(rng);
    timeDerivativePlus[face] = faceDerivatives_ + plusCell * DerivativesSize;
    timeDerivativeMinus[face] = faceDerivatives_ + minusCell * DerivativesSize;

    faceInformation[face].plusSide = sideDist(rng);
    faceInformation[face].minusSide = sideDist(rng);
    // a fault face always addresses the minus side here
    faceInformation[face].faceRelation = 1;
  }

  // the same addresses on the device: at the same offsets into the data there
  if constexpr (isDeviceOn()) {
    auto* drMappingDevice = layer.var<LTS::DRMappingDevice>(Cfg());
    auto* imposedStatePlusDevice = reinterpret_cast<real*>(
        faults.var<DynamicRupture::ImposedStatePlus>(Cfg(), initializer::AllocationPlace::Device));
    auto* fluxSolverPlusDevice = reinterpret_cast<real*>(
        faults.var<DynamicRupture::FluxSolverPlus>(Cfg(), initializer::AllocationPlace::Device));
    real** timeDerivativePlusDevice = faults.var<DynamicRupture::TimeDerivativePlusDevice>(Cfg());
    real** timeDerivativeMinusDevice = faults.var<DynamicRupture::TimeDerivativeMinusDevice>(Cfg());
    for (std::size_t cell = 0; cell < cellCount_; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        const auto& mapping = drMapping[cell][face];
        auto& mappingDevice = drMappingDevice[cell][face];
        mappingDevice.side = mapping.side;
        mappingDevice.faceRelation = mapping.faceRelation;
        mappingDevice.godunov = imposedStatePlusDevice + (mapping.godunov - imposedStatePlus[0]);
        mappingDevice.fluxSolver = fluxSolverPlusDevice + (mapping.fluxSolver - fluxSolverPlus[0]);
      }
    }
    for (std::size_t face = 0; face < faults.size(); ++face) {
      timeDerivativePlusDevice[face] =
          faceDerivativesOnDevice_ + (timeDerivativePlus[face] - faceDerivatives_);
      timeDerivativeMinusDevice[face] =
          faceDerivativesOnDevice_ + (timeDerivativeMinus[face] - faceDerivatives_);
    }
  }

  initFaultState(*factory, frictionLaw);
}

template <typename Cfg>
void SyntheticData<Cfg>::initFaultState(dr::factory::AbstractFactory& factory,
                                        FrictionLaw frictionLaw) {
  using namespace dr::misc::quantity_indices;
  DynamicRupture::Layer& faults = faces();
  constexpr auto PointCount = dr::misc::NumPaddedPoints<Cfg>;

  // the material on both sides
  auto* impAndEta = faults.var<DynamicRupture::ImpAndEta>(Cfg());
  auto* waveSpeedsPlus = faults.var<DynamicRupture::WaveSpeedsPlus>();
  auto* waveSpeedsMinus = faults.var<DynamicRupture::WaveSpeedsMinus>();
  for (std::size_t face = 0; face < faults.size(); ++face) {
    waveSpeedsPlus[face] = FaultMaterial;
    waveSpeedsMinus[face] = FaultMaterial;
    impAndEta[face] = dr::impedancesAndEta<Cfg>(FaultMaterial, FaultMaterial);
  }

  // the initial state, the last stress source of a face, which is in effect from the start (its
  // onset and its rise time stay zero, as those of the nucleations, which have no stress)
  const auto sourceCount = dr::stressSourceCount(*drParameters_);
  const auto initialSource = sourceCount - 1;
  auto* stress = faults.var<DynamicRupture::StressSourceInFaultCS>(Cfg());
  const auto shearStress = shearStressOf(frictionLaw.type);
  for (std::size_t face = 0; face < faults.size(); ++face) {
    for (std::uint32_t point = 0; point < PointCount; ++point) {
      stress[face * sourceCount + initialSource][XX][point] = NormalStress;
      stress[face * sourceCount + initialSource][XY][point] = shearStress;
    }
  }

  // the parameters of the friction law, as the fault files would give them
  const auto initializer = factory.initializer();
  const auto values = faultParametersOf(frictionLaw.type);
  for (const auto& [name, data] : initializer->lawParameters(faults)) {
    const auto value = values.find(name);
    if (value == values.end()) {
      throw std::logic_error("The proxy has no value for the parameter " + name +
                             " of the friction law " + std::to_string(frictionLaw.type) + ".");
    }
    std::fill_n(static_cast<real*>(data), faults.size() * PointCount, value->second);
  }

  // The imposed slip rates take their slip in the coordinate systems of the faces from the mesh
  // otherwise: one meter along the first direction. Their slip began half a rise time before the
  // steps, so that these see the middle of the source time function (the delta of the source time
  // function of FL=35 comes at the start).
  if (frictionLaw.type == 33 || frictionLaw.type == 34 || frictionLaw.type == 35) {
    auto* slipDirection1 = faults.var<LTSImposedSlipRates::ImposedSlipDirection1>(Cfg());
    auto* onsetTime = faults.var<LTSImposedSlipRates::OnsetTime>(Cfg());
    const auto onset = frictionLaw.type == 35 ? 0.0 : -0.5 * values.at("rupture_rise_time");
    for (std::size_t face = 0; face < faults.size(); ++face) {
      std::fill_n(slipDirection1[face], PointCount, 1.0);
      std::fill_n(onsetTime[face], PointCount, onset);
    }
  }

  initializer->deriveInitialState(drStorage_);

  // The initial state of the laws of rate and state holds the faces at rest at their initial slip
  // rate, under the stress it was derived for; with the nucleation of the benchmark on top, the
  // faces rupture instead, as the nucleation patch does in the benchmark.
  const auto nucleationStress = nucleationStressOf(frictionLaw.type);
  for (std::size_t face = 0; face < faults.size(); ++face) {
    for (std::uint32_t point = 0; point < PointCount; ++point) {
      stress[face * sourceCount + initialSource][XY][point] += nucleationStress;
    }
  }
}

template <typename Cfg>
void SyntheticData<Cfg>::initDevice() {
  const bool plasticity = plasticity_ != CellPlasticity::None;
#ifdef ACL_DEVICE
  const auto& device = device::DeviceInstance::instance();
  ltsStorage_.synchronizeTo(initializer::AllocationPlace::Device, device.api().getDefaultStream());
  if (dynamicRupture_) {
    drStorage_.synchronizeTo(initializer::AllocationPlace::Device, device.api().getDefaultStream());
  }
  device.api().syncDefaultStreamWithHost();

  initializer::internal::deriveRequiredScratchpadMemoryForWp(plasticity, ltsStorage_);
  ltsStorage_.allocateScratchPads();
  if (dynamicRupture_) {
    initializer::internal::deriveRequiredScratchpadMemoryForDr(drStorage_);
    drStorage_.allocateScratchPads();
  }
#endif
  initializer::internal::setupRecorders(
      ltsStorage_, drStorage_, plasticity, GravitationalAcceleration);
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg) template class SyntheticData<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

} // namespace seissol::proxy
