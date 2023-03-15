
#ifndef GLOBAL_TIMESTEP_HPP_
#define GLOBAL_TIMESTEP_HPP_
#include <vector>
#include <array>
#include <functional>
#include <Eigen/Dense>

#include "Equations/datastructures.hpp"
#include "Initializer/ParameterDB.h"

namespace seissol::initializer {
    struct GlobalTimestep {
        std::vector<double> elementTimestep;
        double minTimestep;
        double maxTimestep;
    };

    using CellToVertexFunction = std::function<std::array<Eigen::Vector3d, 4>(size_t)>;

    template<typename TQuery, typename TMesh>
    GlobalTimestep computeTimesteps(double cfl, double maximumAllowedTimeStep, const std::string& velocityModel, const TMesh& mesh, size_t cellCount, const CellToVertexFunction& vertices) {
        using Material = seissol::model::MaterialClass;

        TQuery queryGen(mesh);
        std::vector<Material> materials(cellCount);
        seissol::initializers::MaterialParameterDB<Material> parameterDB;
        parameterDB.setMaterialVector(&materials);
        parameterDB.evaluateModel(velocityModel, &queryGen);

        GlobalTimestep timestep;
        timestep.elementTimestep.resize(cellCount);

        for (unsigned cell = 0; cell < cellCount; ++cell) {
            double pWaveVel = materials[cell].getMaxWaveSpeed();

            // Compute insphere radius
            std::array<Eigen::Vector3d, 4> x = vertices(cell);
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
            timestep.elementTimestep[cell] = std::fmin(maximumAllowedTimeStep, cfl * 2.0 * insphere / (pWaveVel * (2 * CONVERGENCE_ORDER - 1)));
        }

        double localMinTimestep = *std::min_element(timestep.elementTimestep.begin(), timestep.elementTimestep.end());
        double localMaxTimestep = *std::max_element(timestep.elementTimestep.begin(), timestep.elementTimestep.end());

        #ifdef USE_MPI
        MPI_Allreduce(&localMinTimestep, &timestep.minTimestep, 1, MPI_DOUBLE, MPI_MIN, seissol::MPI::mpi.comm());
        MPI_Allreduce(&localMaxTimestep, &timestep.maxTimestep, 1, MPI_DOUBLE, MPI_MAX, seissol::MPI::mpi.comm());
        #else
        timestep.minTimestep = localMinTimestep;
        timestep.maxTimestep = localMaxTimestep;
        #endif
        return timestep;
    }
}

#endif
