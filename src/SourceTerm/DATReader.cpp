#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <algorithm>

// Includes to read content of a directory
#include <sys/types.h>
#include <dirent.h>

#include "DATReader.h"



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

	dat->endtime = dat->time[0].back();

}

/*
 * 	Input:
 *  	- position vector3d
 * 		- timestamp
 * 	Output:
 * 		- Sigma_xx value
 * 	Interpolation:
 * 		- Linear in time
 * 		- Inverse Weighted Distance (IWD) in space
 */ 

double seissol::sourceterm::DAT::getSigmaXX(std::vector<double> const& position, double timestamp)
{

	double pos_x = position[0];
	double pos_y = position[1];
	double pos_z = position[2];
	
	
	// IWD with n neighbours
	const int n = 4;

	// From all available receivers hold values needed for IWD with n neighbors
	std::vector<int> rec_idx;			
	std::vector<double> rec_distance;
	std::vector<double> rec_sigma;

	for (int z=0; z < n; ++z) {
		rec_idx.push_back(-1);
		rec_sigma.push_back(0);
		rec_distance.push_back(1e+10);
	}

	for (int i=0; i < pos.size(); ++i) {

		double l2_dist = ((pos[i][0] - pos_x) * (pos[i][0] - pos_x) + (pos[i][1] - pos_y) * (pos[i][1] - pos_y) 
					+ (pos[i][2] - pos_z) * (pos[i][2] - pos_z));

		for (int j=0; j < n; ++j){
			if (l2_dist < rec_distance[j] && !std::count(rec_idx.begin(), rec_idx.end(), i)) {
				rec_distance[j] = l2_dist;
				rec_idx[j] = i;
			}
		}
	}


	// Linear interpolation in time
	// Linear interpolation between two known points (x0, y0) and (x1, y1) for a value x:
    // y = y0 + (x - x0) * (y1 - y0) / (x1 - x0) | Here: x = time , y = sigma_xx
	for (int i=0; i < n; ++i) {

		auto const lower_it = std::lower_bound(time[rec_idx[i]].begin(), time[rec_idx[i]].end(), timestamp);
		auto const upper_it = std::upper_bound(time[rec_idx[i]].begin(), time[rec_idx[i]].end(), timestamp);

		int lower_idx = lower_it - time[rec_idx[i]].begin();
		int upper_idx = upper_it - time[rec_idx[i]].begin();

		if (lower_idx == upper_idx) {
			--lower_idx;
			rec_sigma[i] = sigma_xx[rec_idx[i]][lower_idx] + (timestamp - time[rec_idx[i]][lower_idx]) * (sigma_xx[rec_idx[i]][upper_idx] - sigma_xx[rec_idx[i]][lower_idx])
							/ (time[rec_idx[i]][upper_idx] - time[rec_idx[i]][lower_idx]);
		} else {
			rec_sigma[i] = sigma_xx[rec_idx[i]][timestamp];
		}
	}


	// rec_sigma now stores the linearly interpolated values for all n receivers.
	// Apply IWD with n neighbors (=rec_sigma.size) with power 2

	if (rec_distance[0] == 0) {
		return rec_sigma[rec_idx[0]];
	}

	else {
		double numerator = 0;
		double denominator = 0;

		for (int i=0; i < rec_sigma.size(); ++i){
			double w_i = 1 / (rec_distance[i] * rec_distance[i]);
			numerator += rec_sigma[i] * w_i;
			denominator += w_i;
		}

		return numerator / denominator;
	}

}


// int main( int argc, const char* argv[] )
// {
	
// 	seissol::sourceterm::DAT *dat = new seissol::sourceterm::DAT(); 

// 	double time = 1.5248e+01; // Expected Result: 1.081081636151335e-01 = 0.1081

// 	readDAT( "/Users/philippwendland/Documents/TUM_Master/Semester_4/SeisSol_Results/TRC/output-sin_5e-01Hz_involume",
// 			 dat );


// 	std::vector<double> pos = {5.0, 4.5, 0.1};

// 	double returned_sigma = dat->getSigmaXX(pos, time);

// 	std::cout << returned_sigma << std::endl;


// 	return 0;
// }