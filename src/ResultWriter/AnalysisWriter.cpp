#include "AnalysisWriter.h"

#include <cmath>
#include <string>

#include "SeisSol.h"
#include "Geometry/MeshReader.h"
#include "Physics/InitialField.cpp"

void seissol::writer::AnalysisWriter::printAnalysis(double simulationTime) {
#ifdef USE_MPI
  const auto& mpi = seissol::MPI::mpi;
  const auto& comm = mpi.comm();
#endif // USE_MPI

  const auto initialConditionType = std::string(e_interoperability.getInitialConditionType());
  logInfo(MPI::mpi.rank())
    << "Print analysis for initial conditions" << initialConditionType
    << " at time " << simulationTime;
  
  if (initialConditionType != "Planarwave") {
    return;
  }

  auto* lts = seissol::SeisSol::main.getMemoryManager().getLts();
  auto* ltsLut = e_interoperability.getLtsLut();

  std::vector<Vertex> const& vertices = meshReader->getVertices();
  std::vector<Element> const& elements = meshReader->getElements();
  logInfo(mpi.rank()) << "Elements.size() = " << elements.size();

  using ErrorArray_t = std::array<double, NUMBER_OF_QUANTITIES>;
  auto errL1Local = ErrorArray_t{0.0};
  auto errL2Local = ErrorArray_t{0.0};
  auto errLInfLocal = ErrorArray_t{-1.0};
  auto elemLInfLocal = std::array<unsigned int, NUMBER_OF_QUANTITIES>{0};

  // Initialize quadrature nodes and weights.
  // TODO(Lukas) Increase quadrature order later.
  constexpr auto quadPolyDegree = CONVERGENCE_ORDER+1;
  constexpr auto numQuadPoints = quadPolyDegree * quadPolyDegree * quadPolyDegree;

  double quadraturePoints[numQuadPoints][3];
  double quadratureWeights[numQuadPoints];
  seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, quadPolyDegree);

  // Create basis evaluation objects centered at quad points.
  using BasisSampler = basisFunction::SampledBasisFunctions<real>;
  auto basisVec = std::vector<BasisSampler>{};
  basisVec.reserve(numQuadPoints);
    
  for (int i = 0; i < numQuadPoints; ++i) {
    const auto& curQuadraturePoints = quadraturePoints[i];
    basisVec.push_back(BasisSampler(CONVERGENCE_ORDER,
				    curQuadraturePoints[0],
				    curQuadraturePoints[1],
				    curQuadraturePoints[2]));
  }

  // Note: We iterate over mesh cells by id to avoid
  // cells that are duplicates.
  auto layerMask = initializers::LayerMask(Ghost);
  for (unsigned int meshId = 0; meshId < elements.size(); ++meshId) {
    // Needed to weight the integral.
    const auto volume = MeshTools::volume(elements[meshId], vertices);
    const auto jacobiDet = 6 * volume;

    // Compute global position of quadrature points.
    double const* elementCoords[4];
    double quadraturePointsXyz[numQuadPoints][3];
    for (unsigned v = 0; v < 4; ++v) {
      elementCoords[v] = vertices[elements[meshId].vertices[ v ] ].coords;
    }
    for (unsigned int i = 0; i < numQuadPoints; ++i) {
      seissol::transformations::tetrahedronReferenceToGlobal(elementCoords[0], elementCoords[1], elementCoords[2], elementCoords[3], quadraturePoints[i], quadraturePointsXyz[i]);
    }

    for (size_t i = 0; i < numQuadPoints; ++i) {
      // Evaluate analytical solution at current quad. node
      double analyticalSolution[NUMBER_OF_QUANTITIES];
      const auto x = quadraturePointsXyz[i][0];
      const auto y = quadraturePointsXyz[i][1];
      const auto z = quadraturePointsXyz[i][2];
      initial_field_planarwave(simulationTime, x, y, z, analyticalSolution);

      const auto curWeight = jacobiDet * quadratureWeights[i];
      for (size_t v = 0; v < NUMBER_OF_QUANTITIES; ++v) {
	// Evaluate discrete solution at quad. point
	const auto handle = lts->dofs;
	auto const *coeffsBegin = &ltsLut->lookup(handle, meshId)[v*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
	const auto value = basisVec[i].evalWithCoeffs(coeffsBegin);

	const auto curError = std::abs(value - analyticalSolution[v]);
	errL1Local[v] += curWeight * curError;
	errL2Local[v] += curWeight * curError * curError;

	if (curError > errLInfLocal[v]) {
	  errLInfLocal[v] = curError;
	  elemLInfLocal[v] = meshId;
	}
      }
    }
  }

  for (unsigned int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
    // Find position of element with lowest LInf error.
    VrtxCoords center;
    MeshTools::center(elements[elemLInfLocal[i]],
				    vertices,
				    center);

    /*
    logInfo() << "L1, var[" << i << "] =\t" << errL1Local[i];
    logInfo() << "L2, var[" << i << "] =\t" << std::sqrt(errL2Local[i]);
    logInfo() << "LInf, var[" << i << "] =\t" << errLInfLocal[i]
	      << "\tat [" << center[0] << ",\t" << center[1] << ",\t" << center[2] << "\t]";
    */
  }

  // TODO(Lukas) Print hs, fortran: MESH%MaxSQRTVolume, MESH%MaxCircle
  logInfo() << "Begin reduction";

#ifdef USE_MPI
  // Reduce error over all MPI ranks.
  auto errL1MPI = ErrorArray_t{0.0};
  auto errL2MPI = ErrorArray_t{0.0};
  auto errLInfMPI = ErrorArray_t{0.0};

  MPI_Reduce(errL1Local.data(), errL1MPI.data(), errL1Local.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(errL2Local.data(), errL2MPI.data(), errL2Local.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(errLInfLocal.data(), errLInfMPI.data(), errLInfLocal.size(), MPI_DOUBLE, MPI_MAX, 0, comm);

  if (mpi.rank() == 0) {
    // Log debug output.
    logInfo(mpi.rank()) << "MPI-Error analysis:";
    for (unsigned int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      logInfo(mpi.rank()) << "L1, var[" << i << "] =\t" << errL1MPI[i];
      logInfo(mpi.rank()) << "L2, var[" << i << "] =\t" << std::sqrt(errL2MPI[i]);
      logInfo(mpi.rank()) << "LInf, var[" << i << "] =\t" << errLInfMPI[i];
    }

  }
  MPI_Barrier(comm);
#endif // USE_MPI
}
