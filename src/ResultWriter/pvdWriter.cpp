/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Atanas Atanasov (atanasoa AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Atanas_Atanasov)
 * @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "pvdWriter.h"
#include <iostream>     // std::cout
#include <iomanip>      // std::setw, std::setfill
const std::string PVDWriter::HEADER =
		"<?xml version='1.0' ?>\n "\
		"<VTKFile type='Collection' version='0.1'>\n "\
		"<Collection>\n ";

extern "C"{
#ifdef __xlC__
void create_pvd_writer(long long* ref){
	*ref=(long long) new PVDWriter;
}
void open_pvd_writer(long long* ref){
	((PVDWriter*) *ref)->open();
}
void addtimestep_pvd_writer(long long* ref,double* timestep,int* rank,int* iteration){
	((PVDWriter*) *ref)->addTimestep(*timestep,*rank,*iteration);
}
void close_pvd_writer(long long* ref){
	((PVDWriter*) *ref)->close();
}
void destroy_pvd_writer(long long* ref){
	delete (PVDWriter*) *ref;
}
#else
void create_pvd_writer_(long long* ref, const char* basename){
	*ref=(long long) new PVDWriter(basename);
}

void open_pvd_writer_(long long* ref){
	((PVDWriter*) *ref)->open();
}

void addtimestep_pvd_writer_(long long* ref,double* timestep,int* rank,int* iteration){
	((PVDWriter*) *ref)->addTimestep(*timestep,*rank,*iteration);
}

void close_pvd_writer_(long long* ref){
	((PVDWriter*) *ref)->close();
}

void destroy_pvd_writer_(long long* ref){
	delete (PVDWriter*) *ref;
}
#endif
}

PVDWriter::PVDWriter(const char* basename){
	_baseName = basename;
	_fileName = _baseName + "-fault.pvd";

	// Remove all directories from baseName
	size_t pos = _baseName.find_last_of('/');
	if (pos != std::string::npos)
		_baseName = _baseName.substr(pos+1);
}

PVDWriter::~PVDWriter(){

}

void PVDWriter::open(){

	_out.open( _fileName.c_str() );
	_out<<HEADER<<std::endl;
}

void PVDWriter::addTimestep(double timestep,int rank,int iteration){
	_out<<"<DataSet timestep='"<<timestep<<"' group='' part='"<<rank<<"' file='"<<_baseName<<"-fault-"<<std::setw(5)<<std::setfill('0')<<rank<<"/timestep-"<<std::setw(9)<<std::setfill('0')<<iteration<<".vtu'/>"<<std::endl;
}

void PVDWriter::close(){
	_out<<"</Collection>"<<std::endl;
	_out<<"</VTKFile>"<<std::endl;
	_out.flush();
	_out.close();
}
