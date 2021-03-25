#ifndef SOURCETERM_DATREADER_H_
#define SOURCETERM_DATREADER_H_

#include <vector>
#include <Eigen/Dense>

namespace seissol {
	namespace sourceterm {
		
		struct DAT {
			int q_dim = 6;
			double endtime = 0;
			// IWD with n neighbours
			const int n = 4;

			std::vector<std::vector<double>> pos;
			std::vector<std::vector<double>> time;
			std::vector<std::vector<Eigen::VectorXd>> q;
			
			Eigen::VectorXd getQ(Eigen::Vector3d const& position, double time);

		};


		void readDAT(char const* path, DAT* dat);
	}
}
#endif