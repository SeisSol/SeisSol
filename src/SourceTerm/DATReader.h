#ifndef SOURCETERM_DATREADER_H_
#define SOURCETERM_DATREADER_H_

#include <vector>

namespace seissol {
	namespace sourceterm {
		
		struct DAT {
			std::vector<std::vector<double>> pos;
			std::vector<std::vector<double>> time;
			std::vector<std::vector<double>> sigma_xx;

			double endtime = 0;

			double getSigmaXX(std::vector<double> const& position, double time);

		};


		void readDAT(char const* path, DAT* dat);
	}
}
#endif