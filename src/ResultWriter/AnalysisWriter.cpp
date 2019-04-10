#include "AnalysisWriter.h"

#include <cmath>
#include <string>

#include "SeisSol.h"
#include "Geometry/MeshReader.h"
#include <Physics/InitialField.h>

void seissol::writer::AnalysisWriter::printAnalysis(double simulationTime) {
  const auto& mpi = seissol::MPI::mpi;

  const auto initialConditionType = std::string(e_interoperability.getInitialConditionType());
  logInfo(mpi.rank())
    << "Print analysis for initial conditions" << initialConditionType
    << " at time " << simulationTime;
  
  if (initialConditionType != "Planarwave") {
    return;
  }
  physics::Planarwave iniField;

  auto* lts = seissol::SeisSol::main.getMemoryManager().getLts();
  auto* ltsLut = e_interoperability.getLtsLut();

  std::vector<Vertex> const& vertices = meshReader->getVertices();
  std::vector<Element> const& elements = meshReader->getElements();

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

  std::vector<std::array<double, 3>> quadraturePointsXyz;
  quadraturePointsXyz.resize(numQuadPoints);

  real analyticalSolutionData[numQuadPoints * NUMBER_OF_QUANTITIES] __attribute__((aligned(ALIGNMENT)));
  MatrixView analyticalSolution(analyticalSolutionData, sizeof(analyticalSolutionData)/sizeof(real), &colMjrIndex<numQuadPoints>);

  // Note: We iterate over mesh cells by id to avoid
  // cells that are duplicates.
  for (unsigned int meshId = 0; meshId < elements.size(); ++meshId) {
    // Needed to weight the integral.
    const auto volume = MeshTools::volume(elements[meshId], vertices);
    const auto jacobiDet = 6 * volume;

    // Compute global position of quadrature points.
    double const* elementCoords[4];
    for (unsigned v = 0; v < 4; ++v) {
      elementCoords[v] = vertices[elements[meshId].vertices[ v ] ].coords;
    }
    for (unsigned int i = 0; i < numQuadPoints; ++i) {
      seissol::transformations::tetrahedronReferenceToGlobal(elementCoords[0], elementCoords[1], elementCoords[2], elementCoords[3], quadraturePoints[i], quadraturePointsXyz[i].data());
    }

    // Evaluate analytical solution at quad. nodes
    iniField.evaluate(simulationTime, quadraturePointsXyz, analyticalSolution);

    for (size_t i = 0; i < numQuadPoints; ++i) {
      const auto curWeight = jacobiDet * quadratureWeights[i];
      for (size_t v = 0; v < NUMBER_OF_QUANTITIES; ++v) {
	// Evaluate discrete solution at quad. point
	const auto handle = lts->dofs;
	auto const *coeffsBegin = &ltsLut->lookup(handle, meshId)[v*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
	const auto value = basisVec[i].evalWithCoeffs(coeffsBegin);

	const auto curError = std::abs(value - analyticalSolution(i,v));
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

  }

  // TODO(Lukas) Print hs, fortran: MESH%MaxSQRTVolume, MESH%MaxCircle

#ifdef USE_MPI
  const auto& comm = mpi.comm();

  // Reduce error over all MPI ranks.
  auto errL1MPI = ErrorArray_t{0.0};
  auto errL2MPI = ErrorArray_t{0.0};

  MPI_Reduce(errL1Local.data(), errL1MPI.data(), errL1Local.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(errL2Local.data(), errL2MPI.data(), errL2Local.size(), MPI_DOUBLE, MPI_SUM, 0, comm);


  // Find maximum element and its location.
  auto errLInfSend = std::array<data, errLInfLocal.size()>{};
  auto errLInfRecv = std::array<data, errLInfLocal.size()>{};
  for (size_t i = 0; i < errLInfLocal.size(); ++i) {
    errLInfSend[i] = data{errLInfLocal[i], mpi.rank()};
  }
  MPI_Allreduce(errLInfSend.data(), errLInfRecv.data(),
		errLInfSend.size(),
		MPI_DOUBLE_INT,
		MPI_MAXLOC,
		comm);

  for (unsigned int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
    VrtxCoords centerSend;
    MeshTools::center(elements[elemLInfLocal[i]],
		      vertices,
		      centerSend);


    if (mpi.rank() == errLInfRecv[i].rank &&
	errLInfRecv[i].rank != 0)  {
      MPI_Send(centerSend, 3, MPI_DOUBLE, 0, i, comm);
    }

    if (mpi.rank() == 0) {
      VrtxCoords centerRecv;
      if (errLInfRecv[i].rank == 0) {
	std::copy_n(centerSend, 3, centerRecv);
      } else {
	MPI_Recv(centerRecv, 3, MPI_DOUBLE,  errLInfRecv[i].rank, i, comm, MPI_STATUS_IGNORE);
      }
      logInfo(mpi.rank()) << "L1  , var[" << i << "] =\t" << errL1MPI[i];
      logInfo(mpi.rank()) << "L2  , var[" << i << "] =\t" << std::sqrt(errL2MPI[i]);
      logInfo(mpi.rank()) << "LInf, var[" << i << "] =\t" << errLInfRecv[i].val
			  << "at rank " << errLInfRecv[i].rank
	      << "\tat [" << centerRecv[0] << ",\t" << centerRecv[1] << ",\t" << centerRecv[2] << "\t]";

    }
  }
#else
  for (unsigned int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
    VrtxCoords center;
    MeshTools::center(elements[elemLInfLocal[i]],
		      vertices,
		      center);

    logInfo() << "L1, var[" << i << "] =\t" << errL1Local[i];
    logInfo() << "L2, var[" << i << "] =\t" << std::sqrt(errL2Local[i]);
    logInfo() << "LInf, var[" << i << "] =\t" << errLInfLocal[i]
	      << "\tat [" << center[0] << ",\t" << center[1] << ",\t" << center[2] << "\t]";
  }
#endif // USE_MPI

}

