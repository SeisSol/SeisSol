#ifndef SOURCETERM_DATREADER_H_
#define SOURCETERM_DATREADER_H_

#include <vector>
#include <Eigen/Dense>

namespace seissol {
	namespace sourceterm {
		
		struct DAT {
			std::vector<std::vector<double>> pos;
			std::vector<Eigen::Vector3d> normal;
			std::vector<std::vector<double>> time;
			std::vector<std::vector<double>> pressure_field;
			
			double endtime = 0;
			
			double getPressureField(Eigen::Vector3d const& position, double time);

		};


		int readDAT(char const* path, DAT* dat);
	}
}
#endif