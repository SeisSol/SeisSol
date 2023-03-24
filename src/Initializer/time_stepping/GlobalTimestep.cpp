#include "GlobalTimestep.hpp"

#include <vector>
#include <array>
#include <functional>
#include <Eigen/Dense>

#include "Equations/datastructures.hpp"
#include "Initializer/ParameterDB.h"
#include "Initializer/InputParameters.hpp"

#include "SeisSol.h"

namespace seissol::initializer {
    GlobalTimestep computeTimesteps(double cfl, double maximumAllowedTimeStep, const std::string& velocityModel, const seissol::initializers::C2VArray& ctov) {
        using Material = seissol::model::MaterialClass;

        const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();

        auto* queryGen = seissol::initializers::getBestQueryGenerator(
            seissol::initializer::parameters::modelAnelastic(),
            ssp.model.plasticity,
            seissol::initializer::parameters::modelAnisotropic(),
            seissol::initializer::parameters::modelPoroelastic(),
            ssp.model.useCellHomogenizedMaterial,
            ctov);
        std::vector<Material> materials(ctov.size);
        seissol::initializers::MaterialParameterDB<Material> parameterDB;
        parameterDB.setMaterialVector(&materials);
        parameterDB.evaluateModel(velocityModel, queryGen);

        GlobalTimestep timestep;
        timestep.elementTimeStep.resize(ctov.size);

        for (unsigned cell = 0; cell < ctov.size; ++cell) {
            double pWaveVel = materials[cell].getMaxWaveSpeed();

            // Compute insphere radius
            std::array<Eigen::Vector3d, 4> x = ctov.elementVertices(cell);
            Eigen::Matrix4d A;
            A << x[0](0), x[0](1), x[0](2), 1.0,
                x[1](0), x[1](1), x[1](2), 1.0,
                x[2](0), x[2](1), x[2](2), 1.0,
                x[3](0), x[3](1), x[3](2), 1.0;

            double alpha = A.determinant();
            double Nabc = ((x[1] - x[0]).cross(x[2] - x[0])).norm();
            double Nabd = ((x[1] - x[0]).cross(x[3] - x[0])).norm();
            double Nacd = ((x[2] - x[0]).cross(x[3] - x[0])).norm();
            double Nbcd = ((x[2] - x[1]).cross(x[3] - x[1])).norm();
            double insphere = std::fabs(alpha) / (Nabc + Nabd + Nacd + Nbcd);

            // Compute maximum timestep
            timestep.elementTimeStep[cell] = std::fmin(maximumAllowedTimeStep, cfl * 2.0 * insphere / (pWaveVel * (2 * CONVERGENCE_ORDER - 1)));
        }

        double localMinTimestep = *std::min_element(timestep.elementTimeStep.begin(), timestep.elementTimeStep.end());
        double localMaxTimestep = *std::max_element(timestep.elementTimeStep.begin(), timestep.elementTimeStep.end());

        #ifdef USE_MPI
        MPI_Allreduce(&localMinTimestep, &timestep.globalMinTimeStep, 1, MPI_DOUBLE, MPI_MIN, seissol::MPI::mpi.comm());
        MPI_Allreduce(&localMaxTimestep, &timestep.globalMaxTimeStep, 1, MPI_DOUBLE, MPI_MAX, seissol::MPI::mpi.comm());
        #else
        timestep.globalMinTimeStep = localMinTimestep;
        timestep.globalMaxTimeStep = localMaxTimestep;
        #endif
        return timestep;
    }
}
