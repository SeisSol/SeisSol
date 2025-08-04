// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AnalysisWriter.h"

#include <Alignment.h>
#include <Common/Constants.h>
#include <GeneratedCode/init.h>
#include <GeneratedCode/kernel.h>
#include <GeneratedCode/tensor.h>
#include <Geometry/MeshDefinition.h>
#include <Geometry/MeshTools.h>
#include <Initializer/InitialFieldProjection.h>
#include <Initializer/Parameters/InitializationParameters.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Numerical/Quadrature.h>
#include <Numerical/Transformation.h>
#include <Parallel/MPI.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <mpi.h>
#include <string>
#include <string_view>
#include <utility>
#include <utils/logger.h>
#include <vector>

#include "Geometry/MeshReader.h"
#include "Initializer/PreProcessorMacros.h"
#include "Parallel/OpenMP.h"
#include "Physics/InitialField.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

namespace seissol::writer {

CsvAnalysisWriter::CsvAnalysisWriter(std::string fileName) : fileName(std::move(fileName)) {}

void CsvAnalysisWriter::writeHeader() {
  if (isEnabled) {
    out << "variable,norm,error\n";
  }
}

void CsvAnalysisWriter::addObservation(std::string_view variable,
                                       std::string_view normType,
                                       real error) {
  if (isEnabled) {
    out << variable << "," << normType << "," << error << "\n";
  }
}

void CsvAnalysisWriter::enable() {
  isEnabled = true;
  out.open(fileName);
}

CsvAnalysisWriter::~CsvAnalysisWriter() {
  if (isEnabled) {
    out.close();
    if (!out) {
      logError() << "Error when writing analysis output to file";
    }
  }
}

void AnalysisWriter::printAnalysis(double simulationTime) {
  const auto& mpi = seissol::MPI::mpi;

  const auto initialConditionType = seissolInstance.getSeisSolParameters().initialization.type;
  if (initialConditionType == seissol::initializer::parameters::InitializationType::Zero ||
      initialConditionType == seissol::initializer::parameters::InitializationType::Travelling ||
      initialConditionType ==
          seissol::initializer::parameters::InitializationType::PressureInjection) {
    return;
  }

  logInfo() << "Print analysis for initial conditions" << static_cast<int>(initialConditionType)
            << " at time " << simulationTime;

  const auto& iniFields = seissolInstance.getMemoryManager().getInitialConditions();

  const auto& ltsStorage = seissolInstance.getMemoryManager().getLtsStorage();
  const auto* globalData = seissolInstance.getMemoryManager().getGlobalData().onHost;

  const std::vector<Vertex>& vertices = meshReader->getVertices();
  const std::vector<Element>& elements = meshReader->getElements();

  std::vector<double> data;

  if (initialConditionType == seissol::initializer::parameters::InitializationType::Easi) {
    data = initializer::projectEasiFields(
        {seissolInstance.getSeisSolParameters().initialization.filename},
        simulationTime,
        *meshReader,
        seissolInstance.getSeisSolParameters().initialization.hasTime);
  }

  const auto NumQuantities = 13;

  for (unsigned sim = 0; sim < multisim::NumSimulations; ++sim) {
    logInfo() << "Analysis for simulation" << sim << ": absolute, relative";
    logInfo() << "--------------------------";

    using ErrorArrayT = std::array<double, NumQuantities>;
    using MeshIdArrayT = std::array<unsigned int, NumQuantities>;

    auto errL1Local = ErrorArrayT{0.0};
    auto errL2Local = ErrorArrayT{0.0};
    auto errLInfLocal = ErrorArrayT{-1.0};
    auto elemLInfLocal = MeshIdArrayT{0};
    auto analyticalL1Local = ErrorArrayT{0.0};
    auto analyticalL2Local = ErrorArrayT{0.0};
    auto analyticalLInfLocal = ErrorArrayT{-1.0};

    // also functional with an enabled NVHPC_AVOID_OMP
    const int numThreads = OpenMP::threadCount();
    assert(numThreads > 0);
    // Allocate one array per thread to avoid synchronization.
    auto errsL1Local = std::vector<ErrorArrayT>(numThreads);
    auto errsL2Local = std::vector<ErrorArrayT>(numThreads);
    auto errsLInfLocal = std::vector<ErrorArrayT>(numThreads, {-1});
    auto elemsLInfLocal = std::vector<MeshIdArrayT>(numThreads);
    auto analyticalsL1Local = std::vector<ErrorArrayT>(numThreads);
    auto analyticalsL2Local = std::vector<ErrorArrayT>(numThreads);
    auto analyticalsLInfLocal = std::vector<ErrorArrayT>(numThreads, {-1});

    for (const auto& layer : ltsStorage.leaves(Ghost)) {
      layer.wrap([&](auto cfg) {
        using Cfg = decltype(cfg);
        using real = Real<Cfg>;

      constexpr auto NumQuantities =
          tensor::Q<Cfg>::Shape[sizeof(tensor::Q<Cfg>::Shape) / sizeof(tensor::Q<Cfg>::Shape[0]) - 1];

      // Initialize quadrature nodes and weights.
      // TODO(Lukas) Increase quadrature order later.
      constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
      constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

      double quadraturePoints[NumQuadPoints][3];
      double quadratureWeights[NumQuadPoints];
      seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, QuadPolyDegree);

    // Note: We iterate over mesh cells by id to avoid
    // cells that are duplicates.
    std::vector<std::array<double, 3>> quadraturePointsXyz(NumQuadPoints);

    alignas(Alignment) real numericalSolutionData[tensor::dofsQP<Cfg>::size()];
    alignas(Alignment) real analyticalSolutionData[NumQuadPoints * NumQuantities];

      const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
      const auto* materialData = layer.var<LTS::Material>();
      const auto* dofsData = layer.var<LTS::Dofs>(cfg);
#if defined(_OPENMP) && !NVHPC_AVOID_OMP
      // Note: Adding default(none) leads error when using gcc-8
#pragma omp parallel for shared(elements,                                                          \
                                    vertices,                                                      \
                                    iniFields,                                                     \
                                    quadraturePoints,                                              \
                                    globalData,                                                    \
                                    errsLInfLocal,                                                 \
                                    simulationTime,                                                \
                                    sim,                                                           \
                                    quadratureWeights,                                             \
                                    elemsLInfLocal,                                                \
                                    errsL2Local,                                                   \
                                    errsL1Local,                                                   \
                                    analyticalsL1Local,                                            \
                                    analyticalsL2Local,                                            \
                                    analyticalsLInfLocal)                                          \
    firstprivate(quadraturePointsXyz) private(numericalSolutionData, analyticalSolutionData)
#endif
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        if (secondaryInformation[cell].duplicate > 0) {
          // skip duplicate cells
          continue;
        }
        const auto meshId = secondaryInformation[cell].meshId;
        const int curThreadId = OpenMP::threadId();

        auto numericalSolution = init::dofsQP<Cfg>::view::create(numericalSolutionData);
        auto analyticalSolution = yateto::DenseTensorView<2, real>(analyticalSolutionData,
                                                                   {NumQuadPoints, NumQuantities});

        // Needed to weight the integral.
        const auto volume = MeshTools::volume(elements[meshId], vertices);
        const auto jacobiDet = 6 * volume;

        if (initialConditionType != seissol::initializer::parameters::InitializationType::Easi) {
          // Compute global position of quadrature points.
          const double* elementCoords[Cell::NumVertices];
          for (std::size_t v = 0; v < Cell::NumVertices; ++v) {
            elementCoords[v] = vertices[elements[meshId].vertices[v]].coords;
          }
          for (std::size_t i = 0; i < NumQuadPoints; ++i) {
            seissol::transformations::tetrahedronReferenceToGlobal(elementCoords[0],
                                                                   elementCoords[1],
                                                                   elementCoords[2],
                                                                   elementCoords[3],
                                                                   quadraturePoints[i],
                                                                   quadraturePointsXyz[i].data());
          }

          // Evaluate analytical solution at quad. nodes
          const CellMaterialData& material = materialData[cell];
          iniFields[sim % iniFields.size()]->evaluate(simulationTime,
                                                      quadraturePointsXyz.data(),
                                                      quadraturePointsXyz.size(),
                                                      material,
                                                      analyticalSolution);
        } else {
          for (std::size_t i = 0; i < NumQuadPoints; ++i) {
            for (std::size_t j = 0; j < NumQuantities; ++j) {
              analyticalSolution(i, j) =
                  data.at(meshId * NumQuadPoints * NumQuantities + NumQuantities * i + j);
            }
          }
        }

        auto numSub = seissol::multisim::simtensor(numericalSolution, sim);

        // Evaluate numerical solution at quad. nodes
        kernel::evalAtQP<Cfg> krnl;
        krnl.evalAtQP = globalData->evalAtQPMatrix;
        krnl.dofsQP = numericalSolutionData;
        krnl.Q = dofsData[cell];
        krnl.execute();

        for (size_t i = 0; i < NumQuadPoints; ++i) {
          const auto curWeight = jacobiDet * quadratureWeights[i];
          for (size_t v = 0; v < NumQuantities; ++v) {
            const double curError = std::abs(numSub(i, v) - analyticalSolution(i, v));
            const double curAnalytical = std::abs(analyticalSolution(i, v));

            errsL1Local[curThreadId][v] += curWeight * curError;
            errsL2Local[curThreadId][v] += curWeight * curError * curError;
            analyticalsL1Local[curThreadId][v] += curWeight * curAnalytical;
            analyticalsL2Local[curThreadId][v] += curWeight * curAnalytical * curAnalytical;

            if (curError > errsLInfLocal[curThreadId][v]) {
              errsLInfLocal[curThreadId][v] = curError;
              elemsLInfLocal[curThreadId][v] = meshId;
            }
            analyticalsLInfLocal[curThreadId][v] =
                std::max(curAnalytical, analyticalsLInfLocal[curThreadId][v]);
          }
        }
      }
    });
    }

    for (int i = 0; i < numThreads; ++i) {
      for (unsigned v = 0; v < NumQuantities; ++v) {
        errL1Local[v] += errsL1Local[i][v];
        errL2Local[v] += errsL2Local[i][v];
        analyticalL1Local[v] += analyticalsL1Local[i][v];
        analyticalL2Local[v] += analyticalsL2Local[i][v];
        if (errsLInfLocal[i][v] > errLInfLocal[v]) {
          errLInfLocal[v] = errsLInfLocal[i][v];
          elemLInfLocal[v] = elemsLInfLocal[i][v];
        }
        analyticalLInfLocal[v] = std::max(analyticalsLInfLocal[i][v], analyticalLInfLocal[v]);
      }
    }

    for (std::size_t i = 0; i < NumQuantities; ++i) {
      // Find position of element with lowest LInf error.
      VrtxCoords center;
      MeshTools::center(elements[elemLInfLocal[i]], vertices, center);
    }

    const auto& comm = mpi.comm();

    // Reduce error over all MPI ranks.
    auto errL1MPI = ErrorArrayT{0.0};
    auto errL2MPI = ErrorArrayT{0.0};
    auto analyticalL1MPI = ErrorArrayT{0.0};
    auto analyticalL2MPI = ErrorArrayT{0.0};

    MPI_Reduce(errL1Local.data(), errL1MPI.data(), errL1Local.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(errL2Local.data(), errL2MPI.data(), errL2Local.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(analyticalL1Local.data(),
               analyticalL1MPI.data(),
               analyticalL1Local.size(),
               MPI_DOUBLE,
               MPI_SUM,
               0,
               comm);
    MPI_Reduce(analyticalL2Local.data(),
               analyticalL2MPI.data(),
               analyticalL2Local.size(),
               MPI_DOUBLE,
               MPI_SUM,
               0,
               comm);

    // Find maximum element and its location.
    auto errLInfSend = std::array<Data, errLInfLocal.size()>{};
    auto errLInfRecv = std::array<Data, errLInfLocal.size()>{};
    for (size_t i = 0; i < errLInfLocal.size(); ++i) {
      errLInfSend[i] = Data{errLInfLocal[i], mpi.rank()};
    }
    MPI_Allreduce(errLInfSend.data(),
                  errLInfRecv.data(),
                  errLInfSend.size(),
                  MPI_DOUBLE_INT,
                  MPI_MAXLOC,
                  comm);

    auto analyticalLInfMPI = ErrorArrayT{0.0};
    MPI_Reduce(analyticalLInfLocal.data(),
               analyticalLInfMPI.data(),
               analyticalLInfLocal.size(),
               MPI_DOUBLE,
               MPI_MAX,
               0,
               comm);

    auto csvWriter = CsvAnalysisWriter(fileName);

    if (mpi.rank() == 0) {
      csvWriter.enable();
      csvWriter.writeHeader();
    }

    for (std::size_t i = 0; i < NumQuantities; ++i) {
      VrtxCoords centerSend{};
      MeshTools::center(elements[elemLInfLocal[i]], vertices, centerSend);

      if (mpi.rank() == errLInfRecv[i].rank && errLInfRecv[i].rank != 0) {
        MPI_Send(centerSend, 3, MPI_DOUBLE, 0, i, comm);
      }

      if (mpi.rank() == 0) {
        VrtxCoords centerRecv{};
        if (errLInfRecv[i].rank == 0) {
          std::copy_n(centerSend, 3, centerRecv);
        } else {
          MPI_Recv(centerRecv, 3, MPI_DOUBLE, errLInfRecv[i].rank, i, comm, MPI_STATUS_IGNORE);
        }

        const auto errL1 = errL1MPI[i];
        const auto errL2 = std::sqrt(errL2MPI[i]);
        const auto errLInf = errLInfRecv[i].val;
        const auto errL1Rel = errL1 / analyticalL1MPI[i];
        const auto errL2Rel = std::sqrt(errL2MPI[i] / analyticalL2MPI[i]);
        const auto errLInfRel = errLInf / analyticalLInfMPI[i];
        logInfo() << "L1  , var[" << i << "] =\t" << errL1 << "\t" << errL1Rel;
        logInfo() << "L2  , var[" << i << "] =\t" << errL2 << "\t" << errL2Rel;
        logInfo() << "LInf, var[" << i << "] =\t" << errLInf << "\t" << errLInfRel << "at rank "
                  << errLInfRecv[i].rank << "\tat [" << centerRecv[0] << ",\t" << centerRecv[1]
                  << ",\t" << centerRecv[2] << "\t]";
        csvWriter.addObservation(std::to_string(i), "L1", errL1);
        csvWriter.addObservation(std::to_string(i), "L2", errL2);
        csvWriter.addObservation(std::to_string(i), "LInf", errLInf);
        csvWriter.addObservation(std::to_string(i), "L1_rel", errL1Rel);
        csvWriter.addObservation(std::to_string(i), "L2_rel", errL2Rel);
        csvWriter.addObservation(std::to_string(i), "LInf_rel", errLInfRel);
      }
    }
  }
}
} // namespace seissol::writer
