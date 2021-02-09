#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>

// Includes to read content of a directory
#include <sys/types.h>
#include <dirent.h>

#include "DATReader.h"

// std::vector<std::vector<double>> seissol::sourceterm::DAT::pos;
// std::vector<std::vector<double>> seissol::sourceterm::DAT::time;
// std::vector<std::vector<double>> seissol::sourceterm::DAT::sigma_xx;



void seissol::sourceterm::readDAT(char const* path, DAT* dat)
{
	// Receiver data files always end on '.dat'
	const std::string ending = ".dat";

	// Read content inside path 
	std::vector<std::string> v;
	std::string path_str = path;

	DIR* dirp = opendir(path);
    struct dirent *dp;
    while ((dp = readdir(dirp)) != NULL) {
		std::string current_file_name = dp->d_name;

		if (current_file_name.length() > ending.length() && current_file_name.compare (current_file_name.length() - ending.length(), ending.length(), ending) == 0) 
	    	v.push_back(current_file_name);
    }
    closedir(dirp);


	for (const auto& filename : v) {
		// std::cout << filename << std::endl;

		std::string file_path = path_str + "/" + filename;


		// Save values from current file
		std::vector<double> current_pos;
		std::vector<double> current_time;
		std::vector<double> current_sigmaxx;


		std::ifstream dat_file(file_path);

		assert(dat_file.is_open());
		

		std::string line;

		while (!dat_file.eof())
		{
			getline(dat_file, line, '\n');

			// Lines with position information start with '#'
			if (line.find('#') == 0) {
				std::smatch match;
				// Matches position information of .dat files
				if ( std::regex_search(line, match, std::regex("[-|\\d][0-9.]+e[-+][0-9]+"))) {
					std::string next_pos = match[0];
					current_pos.push_back(std::stod(next_pos));
				}
			// All other important lines start with ' ' (white space)
			} else if ((line.find(' ') == 0)) {

				std::istringstream iss(line);
				std::vector<std::string> tmp;

				for (std::string ex_line; iss >> ex_line; ) {
					tmp.push_back(ex_line);

				}

				current_time.push_back(std::stod(tmp.at(0)));
				current_sigmaxx.push_back(std::stod(tmp.at(1)));
			}
		}

		// Save current data 
		dat->pos.push_back(current_pos);
		dat->time.push_back(current_time);
		dat->sigma_xx.push_back(current_sigmaxx);

	}

	endtime = dat->pos[0].back();


	assert(dat->sigma_xx.size() != 0);



	std::cout << "Reading data successful\n";

}

/*
 *	Return sigma_xx at the given timestamp using linear interpolation.
 *	For now: nearest neighbor
 */ 

double seissol::sourceterm::DAT::getSigmaXX(std::vector<double> const& position, double timestamp)
{

	double pos_x = position[0];
	double pos_y = position[1];
	double pos_z = position[2];
	

	double min_dist = 1e+10;
	int min_idx = -1;

	for (int i = 0; i < pos.size(); ++i) {
		double l2_dist = ((pos[i][0] - pos_x) * (pos[i][0] - pos_x) + (pos[i][1] - pos_y) * (pos[i][1] - pos_y) 
					+ (pos[i][2] - pos_z) * (pos[i][2] - pos_z));
		
		if (l2_dist < min_dist ) {
			min_dist = l2_dist;
			min_idx = i;
		}
	}


	auto const it = std::lower_bound(time[min_idx].begin(), time[min_idx].end(), timestamp);
	if (it == time[min_idx].end()) {return -1;}

	int idx = it - time[min_idx].begin();



	// std::cout << "Extracted index: " << idx << std::endl;
	// std::cout << "Extracted time: " << time[min_idx][idx] << std::endl;
	// std::cout << "Extracted sigma: " << sigma_xx[min_idx][idx] << std::endl;


	// std::cout << sigma_xx[min_idx][idx] << '\n';

	return sigma_xx[min_idx][idx];

}


// int main( int argc, const char* argv[] )
// {
	
// 	seissol::sourceterm::DAT dat; 

// 	double time = 1.039000000000014e+01;

// 	readDAT( "/Users/philippwendland/Documents/TUM_Master/Semester_4/SeisSol_Results/cube_forward/output-sin4",
// 			 dat );


// 	std::vector<double> pos = {-4.5, 4.5, 0.1};

// 	double returned_sigma = dat.getSigmaXX(pos, time);

// 	std::cout << returned_sigma << std::endl;

// 	return 0;
// }