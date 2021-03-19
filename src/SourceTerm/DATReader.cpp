#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <algorithm>
#include <Eigen/Dense>
#include <math.h>
// Includes to read content of a directory
#include <sys/types.h>
#include <dirent.h>

#include "DATReader.h"



/*
* Given a path, readDAT reads all recevier files and a .ndat file containing the normal vector for each receiver 
* on the surface.
* The scalar pressure field for each receiver is saved. It is obtained by extracting the stress tensor at a given
* position and multiplying it by the normal vector and finally computing the length of the resulting vector.
*
*/
int seissol::sourceterm::readDAT(char const* path, DAT* dat)
{
	// Receiver data files always end on '.dat'
	const std::string ending = ".dat";
	
	// the corresponding normal vector for each receiver is saved in a .ndat file
	const std::string normal_vec_ending = ".ndat";
	std::string normal_file_name = "null.ndat";

	// Read content inside path 
	std::vector<std::string> dat_file_list;
	std::string path_str = path;

	DIR* dirp = opendir(path);
    struct dirent *dp;
    while ((dp = readdir(dirp)) != NULL) {
		std::string current_file_name = dp->d_name;
		if (current_file_name.length() > normal_vec_ending.length() && 
			current_file_name.compare (current_file_name.length() - normal_vec_ending.length(), normal_vec_ending.length(), normal_vec_ending) == 0) 
				normal_file_name = current_file_name;

		if (current_file_name.length() > ending.length() && 
			current_file_name.compare (current_file_name.length() - ending.length(), ending.length(), ending) == 0) 
	    		dat_file_list.push_back(current_file_name);
    
	}
    closedir(dirp);


	// Save normal vectors as a 3D vector
	// Entry 0 in tmp_normal_vectors stores the normal vector of receiver 1, etc
	std::vector<Eigen::Vector3d> tmp_normal_vectors;
	int skip_rows = 0;

	std::ifstream normal_file(path_str + "/" + normal_file_name);

	if (normal_file.is_open()) {
		std::string normal_line;
		int line_count = 0;
		while(!normal_file.eof()){
			getline(normal_file, normal_line, '\n');
			
			std::istringstream iss(normal_line);
			std::vector<std::string> tmp;
			for (std::string ex_line; iss >> ex_line; ) {
				tmp.push_back(ex_line);
				}
			if (line_count == 0) {
				skip_rows = std::stod(tmp.at(0));
			} else {
				Eigen::Vector3d eigen_tmp(std::stod(tmp.at(0)), std::stod(tmp.at(1)), std::stod(tmp.at(2)));
				tmp_normal_vectors.push_back(eigen_tmp);
			}

			++line_count;
		}
	} else {
		std::cerr << "Could not process	.ndat file \n";
		return 1;
	}

	for (const auto& filename : dat_file_list) {

		std::string file_path = path_str + "/" + filename;
		
		// rec_nr will hold the number of the current receiver file
		// needed in order to extract the appropriate normal vector
		std::smatch file_name_match;
		std::regex_search(filename, file_name_match, std::regex("-\\d\\d\\d\\d\\d-"));

		std::string file_number = file_name_match[0];
		file_number = file_number.substr(1, file_number.size() - 2);
		int rec_nr = std::stoi(file_number);

		// Temporarly save values from current file
		std::vector<double> current_pos;
		Eigen::Vector3d current_normal;
		std::vector<double> current_time;
		std::vector<double> current_pressure_field;

		current_normal = tmp_normal_vectors.at(rec_nr - skip_rows - 1);

		std::ifstream dat_file(file_path);
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
				
				Eigen::MatrixXd stress_tensor(3, 3);
				
				// Diagonal elements
				stress_tensor(0, 0) = std::stod(tmp.at(1));
				stress_tensor(1, 1) = std::stod(tmp.at(2));
				stress_tensor(2, 2) = std::stod(tmp.at(3));
				// off-diagonal elements
				stress_tensor(0, 1) = std::stod(tmp.at(4));
				stress_tensor(0, 2) = std::stod(tmp.at(6));
				stress_tensor(1, 0) = stress_tensor(0, 1);
				stress_tensor(2, 0) = stress_tensor(0, 2);
				stress_tensor(1, 2) = std::stod(tmp.at(5));
				stress_tensor(2, 1) = stress_tensor(1, 2);
				
				Eigen::Vector3d normal_stress = stress_tensor * current_normal;

				double p_res = sqrt( normal_stress.dot(normal_stress) );

				current_pressure_field.push_back(p_res);
			}
		}

		// Save current data 
		dat->pos.push_back(current_pos);
		dat->normal.push_back(current_normal);
		dat->time.push_back(current_time);
		dat->pressure_field.push_back(current_pressure_field);

	}

	dat->endtime = dat->time[0].back();

	return 0;

}

/*
 * 	Input:
 *  	- position vector3d
 * 		- timestamp
 * 	Output:
 * 		- Scalar value of pressure field
 * 	Interpolation:
 * 		- Linear in time
 * 		- Inverse Weighted Distance (IWD) in space
 */ 

double seissol::sourceterm::DAT::getPressureField(Eigen::Vector3d const& position, double timestamp)
{
	// timestamp < 0 --> past recording time of forward receivers.
	if (timestamp < 0) {
		return 0;
	}

	double pos_x = position(0);
	double pos_y = position(1);
	double pos_z = position(2);
	
	
	// IWD with n neighbours
	const int n = 4;

	// From all available receivers hold values needed for IWD with n neighbors
	std::vector<int> rec_idx;			
	std::vector<double> rec_distance;
	std::vector<double> rec_sigma;

	for (int z=0; z < n; ++z) {
		rec_idx.push_back(-1);
		rec_sigma.push_back(0);
		rec_distance.push_back(1e+20);
	}

	// Use L2-norm to find the n receivers closest to the current position
	for (int i=0; i < pos.size(); ++i) {
		double l2_dist = ((pos[i][0] - pos_x) * (pos[i][0] - pos_x) + (pos[i][1] - pos_y) * (pos[i][1] - pos_y) 
					+ (pos[i][2] - pos_z) * (pos[i][2] - pos_z));

		int idx_tmp = std::max_element(rec_distance.begin(), rec_distance.end()) - rec_distance.begin();

		if (rec_distance[idx_tmp] > l2_dist) {
			rec_distance[idx_tmp] = l2_dist;
			rec_idx[idx_tmp] = i;
		}
	}


	// Linear interpolation in time
	// Linear interpolation between two known points (x0, y0) and (x1, y1) for a value x:
    // y = y0 + (x - x0) * (y1 - y0) / (x1 - x0) | Here: x = time , y = pressure_field
	for (int i=0; i < n; ++i) {
		auto const lower_it = std::lower_bound(time[rec_idx[i]].begin(), time[rec_idx[i]].end(), timestamp);
		auto const upper_it = std::upper_bound(time[rec_idx[i]].begin(), time[rec_idx[i]].end(), timestamp);

		int lower_idx = lower_it - time[rec_idx[i]].begin();
		int upper_idx = upper_it - time[rec_idx[i]].begin();

		if (lower_idx == upper_idx) {
			if (lower_idx == 0) {
				upper_idx = 1;
			} else {
				--lower_idx;
			}
			rec_sigma[i] = pressure_field[rec_idx[i]][lower_idx] + (timestamp - time[rec_idx[i]][lower_idx]) * (pressure_field[rec_idx[i]][upper_idx] - pressure_field[rec_idx[i]][lower_idx])
							/ (time[rec_idx[i]][upper_idx] - time[rec_idx[i]][lower_idx]);
		} else {
			rec_sigma[i] = pressure_field[rec_idx[i]][timestamp];
		}
	}

	// rec_sigma now stores the linearly interpolated values for all n receivers.
	// Apply IWD with n neighbors (=rec_sigma.size) with power 2
  	
	int idx_tmp = std::min_element(rec_distance.begin(), rec_distance.end()) - rec_distance.begin();

	if (rec_distance[idx_tmp] == 0.0) {
		return rec_sigma[idx_tmp];
	} else {
		double numerator = 0;
		double denominator = 0;
		for (int i=0; i < n; ++i){
			double w_i = 1 / (rec_distance[i] * rec_distance[i]);
			numerator += rec_sigma[i] * w_i;
			denominator += w_i;
		}
		
		return (numerator / denominator);	
	}
}


// int main( int argc, const char* argv[] )
// {
	
// 	seissol::sourceterm::DAT *dat = new seissol::sourceterm::DAT(); 

// 	double time = 2.47+01; // Expected Result: 1.081081636151335e-01 = 0.1081

// 	// readDAT( "/Users/philippwendland/Documents/TUM_Master/Semester_4/SeisSol_Results/cube_forward/output_7",
// 	// 		 dat );

	
// 	readDAT( "/Users/philippwendland/Documents/TUM_Master/Semester_4/SeisSol_Results/point-source/output-fwd-500kmesh",
// 			 dat );

// 	Eigen::Vector3d pos(-5, 4.5, -0.1);

// 	double returned_p = dat->getPressureField(pos, time);

// 	std::cout << returned_p << std::endl;


// 	return 0;
// }
