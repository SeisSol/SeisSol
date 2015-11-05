/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Atanas Atanasov (atanasoa AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Atanas_Atanasov)
 * @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
 * @author Thomas Ulrich (ulrich AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/ulrich)
 * see for information about binary format http://www.vtk.org/Wiki/VTK_XML_Formats
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
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

#include "vtkWriter.h"
#include <sys/stat.h>   // mkdir()
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setw, std::setfill
#include <fstream>
#include <climits>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
using namespace std;

/*Alice TODO remove hard coding of BigEndian */
const std::string vtkWriter::HEADER = "<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian' header_type='UInt32'>\n" \
"<UnstructuredGrid>\n" \

;

extern "C" {
#ifdef __xlC__
void create_vtk_writer(long long* ref,int* rank,int* iteration){
	*ref=(long long)new vtkWriter(*rank,*iteration);
}

void destroy_vtk_writer(long long* ref){
	delete	(vtkWriter*)*ref;
}

void insert_vertex_vtk_writer(long long* ref,double* x,double* y,double* z){
	((vtkWriter*)*ref)->insert_vertex(*x,*y,*z);
}
void write_vertices_vtk_writer(long long* ref){
	((vtkWriter*)*ref)->write_vertices();
}

void start_cell_data_vtk_writer(long long* ref, int* varid){
	((vtkWriter*)*ref)->start_cell_data(*varid);
}

void end_cell_data_vtk_writer(long long* ref){
	((vtkWriter*)*ref)->end_cell_data();
}

void plot_cell_data_vtk_writer(long long* ref,double* value){
	((vtkWriter*)*ref)->plot_cell_data(*value);
}

void plot_cells_vtk_writer(long long* ref){
	((vtkWriter*)*ref)->plot_cells();
}

void open_vtk_writer(long long* ref){
	((vtkWriter*)*ref)->open();
}

void close_vtk_writer(long long* ref){
	((vtkWriter*)*ref)->close();
}
#else
void create_vtk_writer_(long long* ref,int* rank,int* iteration, const char* basename, int* binaryoutput){
	*ref=(long long)new vtkWriter(*rank,*iteration, basename, *binaryoutput);
}

void destroy_vtk_writer_(long long* ref){
	delete	(vtkWriter*)*ref;
}

void insert_vertex_vtk_writer_(long long* ref,double* x,double* y,double* z){
	((vtkWriter*)*ref)->insert_vertex(*x,*y,*z);
}

void write_vertices_vtk_writer_(long long* ref){
	((vtkWriter*)*ref)->write_vertices();
}

void start_cell_data_vtk_writer_(long long* ref, int* varid){
	((vtkWriter*)*ref)->start_cell_data(*varid);
}

void end_cell_data_vtk_writer_(long long* ref){
	((vtkWriter*)*ref)->end_cell_data();
}

void plot_cell_data_vtk_writer_(long long* ref,double* value){
	((vtkWriter*)*ref)->plot_cell_data(*value);
}

void plot_cells_vtk_writer_(long long* ref){
	((vtkWriter*)*ref)->plot_cells();
}

void open_vtk_writer_(long long* ref){
	((vtkWriter*)*ref)->open();
}

void close_vtk_writer_(long long* ref){
	((vtkWriter*)*ref)->close();
}
#endif
}

vtkWriter::vtkWriter(int rank,int iteration, const char* basename, int binaryoutput)
	: _dirCreated(false) {
	if (!_dirCreated) {
		// Create directory for this rank
		std::stringstream dirname;
		dirname << basename<<"-fault-"<<std::setw(5)<<std::setfill('0')<<rank;
		mkdir(dirname.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);

		_dirCreated = true;
	}

	std::stringstream filename;
	filename<<basename<<"-fault-"<<std::setw(5)<<std::setfill('0')<<rank<<"/timestep-"<<std::setw(9)<<std::setfill('0')<<iteration<<".vtu";
	_fileName=filename.str();
	_number_of_cells=0;
	_vertex_counter=0;
        m_offset=0;
        m_iCurrent=0;
        m_var_id=0;
        //if m_bBinary=true binary file format considered, else ascii
        //if m_bFloat64=true, double written, else float
        if (binaryoutput==0) {
           m_bBinary=false;
        } else if (binaryoutput==1) {
           m_bBinary=true;
           m_bFloat64=false;
        } else {
           m_bBinary=true;
           m_bFloat64=true;
        }

        m_nvar = 15;
        bMagnitudeWritten = new bool [m_nvar];

	for(int i=0;i<m_nvar;i++)
           bMagnitudeWritten[i]=false;
}

vtkWriter::~vtkWriter(){
        delete[] bMagnitudeWritten;
        if (m_bBinary) {
           if (m_bFloat64) 
              delete[] dArraySca;
           else
              delete[] fArraySca;
           }
}

void vtkWriter::insert_vertex(double x,double y,double z){
	int id=0;
	_map.get(x,y,z,id);
	if(id>=0){
		_vertex_counter++;

		_v_x.push_back(x);
		_v_y.push_back(y);
		_v_z.push_back(z);
	}
	_number_of_cells++;

}
void vtkWriter::write_vertices(){
	_out<<"<Points>"<<std::endl;
        if (m_bBinary) {
           if (m_bFloat64) {
	      _out<<"<DataArray type='Float64' NumberOfComponents='3' format='append' offset='"<<m_offset<<"'/>"<<std::endl;
              m_offset = m_offset + sizeof(int) + 3*_vertex_counter*sizeof(double);
              dArraySca = new double [m_nvar*_number_of_cells/3];
           } else {
	      _out<<"<DataArray type='Float32' NumberOfComponents='3' format='append' offset='"<<m_offset<<"'/>"<<std::endl;
              m_offset = m_offset + sizeof(int) + 3*_vertex_counter*sizeof(float);
              fArraySca = new float [m_nvar*_number_of_cells/3];
           }
        } else {
           _out<<"<DataArray type='Float64' NumberOfComponents='3' Format='ascii'>"<<std::endl;
	   for(int i=0;i<_vertex_counter;i++)
	      _out<< _v_x[i]<<" "<<_v_y[i]<<" "<<_v_z[i]<<std::endl;
	   _out<<"</DataArray>"<<std::endl;
        }
	_out<<"</Points>"<<std::endl;
}
void vtkWriter::start_cell_data(int var_id){
	std::string labels [m_nvar];
	labels[0]="SRs";
	labels[1]="SRd";
	labels[2]="T_s";
	labels[3]="T_d";
	labels[4]="P_n";
	labels[5]="u_n";
	labels[6]="Mud";
	labels[7]="StV";
	labels[8]="Ts0";
	labels[9]="Td0";
	labels[10]="Pn0";
	labels[11]="Sls";
	labels[12]="Sld";
	labels[13]="Vr";
	labels[14]="ASl";
        if (m_bBinary) {
           bMagnitudeWritten[var_id-1]=true;
           m_iCurrent=0;
           m_var_id=var_id-1;
           if (m_bFloat64) {
	      _out<<"<DataArray type='Float64' Name='"<<labels[var_id-1].c_str()<<"' format='append'  offset='"<<m_offset<<"'/>"<<std::endl;
              m_offset = m_offset + (sizeof(int) + _number_of_cells*sizeof(double)/3);
           } else {
	      _out<<"<DataArray type='Float32' Name='"<<labels[var_id-1].c_str()<<"' format='append'  offset='"<<m_offset<<"'/>"<<std::endl;
              m_offset = m_offset + (sizeof(int) + _number_of_cells*sizeof(float)/3);
          }
        } else {
          _out<<"<DataArray type='Float64' Name='"<<labels[var_id-1].c_str()<<"' Format='ascii'>"<<std::endl;
        }
}
void vtkWriter::end_cell_data(){
        if (!m_bBinary) 
	   _out<<"</DataArray>"<<std::endl;
}
void vtkWriter::plot_cell_data(double value){
        if (m_bBinary) {
           if (m_bFloat64)
              dArraySca[m_var_id*_number_of_cells/3 + m_iCurrent]=value;
           else
              fArraySca[m_var_id*_number_of_cells/3 + m_iCurrent]=(float) value;
           m_iCurrent++;  
        } else {
	   _out<<value<<std::endl;
        }
}
void vtkWriter::plot_cells(){
	_out<<"<Cells>"<<std::endl;
        if (m_bBinary) {
	   _out<<"<DataArray type='Int32' Name='connectivity' format='append'  offset='"<<m_offset<<"'/>"<<std::endl;
           m_offset = m_offset + sizeof(int) +_number_of_cells*sizeof(int);
        } else {
           _out<<"<DataArray type='Int32' Name='connectivity' Format='ascii'>"<<std::endl;
        
	   int id_v0=0;
	   int id_v1=0;
	   int id_v2=0;
	    /* Alice: for very large meshes (large number of vertices), this might have to be changed to int64, long int */
	   for(int i=0;i<_number_of_cells/3;i++){
	   	_map.get_vid(i*3,id_v0);
	   	_map.get_vid(i*3+1,id_v1);
	   	_map.get_vid(i*3+2,id_v2);
	   	   _out<<id_v0<<" "<<id_v1<<" "<<id_v2<<std::endl;
	   }
	   _out<<"</DataArray>"<<std::endl;
	}
        if (m_bBinary) {
	   _out<<"<DataArray type='Int32' Name='offsets' format='append'  offset='"<<m_offset<<"'/>"<<std::endl;
           m_offset = m_offset + sizeof(int) +_number_of_cells*sizeof(int)/3;
        } else {
           _out<<"<DataArray type='Int32' Name='offsets' Format='ascii'>"<<std::endl;
	   for(int i=0;i<_number_of_cells/3;i++)
	   	_out<<(i+1)*3<<std::endl;
	   _out<<"</DataArray>"<<std::endl;
        }

        if (m_bBinary) {
	   _out<<"<DataArray type='Int32' Name='types' format='append'  offset='"<<m_offset<<"'/>"<<std::endl;
           m_offset = m_offset + sizeof(int) +_number_of_cells*sizeof(int)/3;
        } else {
           _out<<"<DataArray type='Int32' Name='types' Format='ascii'>"<<std::endl;
	   for(int i=0;i<_number_of_cells/3;i++)
	   	_out<<"5"<<std::endl;
	   _out<<"</DataArray>"<<std::endl;
        }
	_out<<"</Cells>"<<std::endl;
	_out<<"<CellData Scalars='scalars'>"<<std::endl;
}
void vtkWriter::open(){
	_out.open( _fileName.c_str(), ios::binary);
	_out<<HEADER<<std::endl;
	_out<<"<Piece NumberOfPoints='"<<_vertex_counter<<"' NumberOfCells='"<<_number_of_cells/3<<"'>"<<std::endl;

}
void vtkWriter::close(){
	_out<<"</CellData>"<<std::endl;
	_out<<"</Piece>"<<std::endl;
	_out<<"</UnstructuredGrid>"<<std::endl;
	if (m_bBinary) {
        // Once supported by Paraview, larger ByteCountType (example long) can be set here
        // 'header_type' have also to be modified
        typedef unsigned int tByteCount;
        unsigned long long MaxInt= UINT_MAX;
 
        tByteCount iByteCount;
        unsigned long long MaxByte;        

        if (m_bFloat64) {
           MaxByte = max(max(3*_vertex_counter*sizeof(double), _number_of_cells*sizeof(int)), _number_of_cells*sizeof(double)/3);
        } else {
           MaxByte = max(max(3*_vertex_counter*sizeof(float), _number_of_cells*sizeof(int)), _number_of_cells*sizeof(float)/3);
        }
        if (MaxByte>MaxInt) { 
           std::cerr<<"to many data for the Paraview binary format "<<MaxByte<<" Bytes, greater than MaxInt "<<MaxInt<<"\n";
           exit(-1);
        }        

	_out<<"<AppendedData encoding='raw'>"<<std::endl<<"_";

        // coordinates
        if (m_bFloat64) {
           iByteCount=3*_vertex_counter*sizeof(double);
           _out.write((const char*) &iByteCount, sizeof(tByteCount));
           double* dArray=new double [3];
           for(int i=0;i<_vertex_counter;i++) {
                   dArray[0]=_v_x[i];
                   dArray[1]=_v_y[i];
                   dArray[2]=_v_z[i];
                   _out.write((const char*) dArray, 3*sizeof(double));
           }        
        delete[] dArray;
        } else {
            iByteCount=3*_vertex_counter*sizeof(float);
           _out.write((const char*) &iByteCount, sizeof(tByteCount));
           float* fArray=new float [3];
           for(int i=0;i<_vertex_counter;i++) {
                   fArray[0]=(float) _v_x[i];
                   fArray[1]=(float) _v_y[i];
                   fArray[2]=(float) _v_z[i];
                   _out.write((const char*) fArray, 3*sizeof(float));
           }        
        delete[] fArray;
        }
     
        // connectivity
        iByteCount = _number_of_cells*sizeof(int);
        _out.write((const char*) &iByteCount, sizeof(tByteCount));
        int* iArray=new int [3];
	for(int i=0;i<_number_of_cells/3;i++){
		_map.get_vid(i*3,  iArray[0]);
		_map.get_vid(i*3+1,iArray[1]);
		_map.get_vid(i*3+2,iArray[2]);
        	_out.write((const char*) iArray, 3*sizeof(int));
	}
        delete[] iArray;
        // offset
        iByteCount = _number_of_cells*sizeof(int)/3;
        _out.write((const char*) &iByteCount, sizeof(tByteCount));
	for(int i=0;i<_number_of_cells/3;i++) {
                int ival =(i+1)*3;
        	_out.write((const char*) &ival, sizeof(int));
        }

        // type
        iByteCount = _number_of_cells*sizeof(int)/3;
        _out.write((const char*) &iByteCount, sizeof(tByteCount));
	for(int i=0;i<_number_of_cells/3;i++) {
                int ival =5;
        	_out.write((const char*) &ival, sizeof(int));
        }
        // SR, SD, etc
        if (m_bFloat64) {
           iByteCount = _number_of_cells*sizeof(double)/3;
	   for(int i=0;i<m_nvar;i++) {
              if (bMagnitudeWritten[i]==true) {
                 _out.write((const char*) &iByteCount, sizeof(tByteCount));
                 double *dSubArray = dArraySca + i*_number_of_cells/3;
                 _out.write((const char*) dSubArray, _number_of_cells/3*sizeof(double));
              }
           }
        } else {
           iByteCount = _number_of_cells*sizeof(float)/3;
	   for(int i=0;i<m_nvar;i++) {
              if (bMagnitudeWritten[i]==true) {
                 _out.write((const char*) &iByteCount, sizeof(tByteCount));
                 float *fSubArray = fArraySca + i*_number_of_cells/3;
                 _out.write((const char*) fSubArray, _number_of_cells/3*sizeof(float));
              }
           }
        }
	_out<<"</AppendedData>"<<std::endl;
        }
	_out<<"</VTKFile>"<<std::endl;
	_out.flush();
	_out.close();
}

/*Alice (comments, to be edited/removed as proceeding)
 *  towards output as appended binary raw, not encoded 		   
 *  encoded binary would require writing a 64-bit encoder (or depending on vtk library) 
 * Structure:
 *  0 ) each XML ataArray tag needs a byte offset entry 
 *  --> The offset variable must count both the binary integer used to represent the length of the data + the actual data
 *  --> the offset has to be summed up after each field
 *  --> the binary integer of array length has to be specified as 32 or 64 
 *  --> Thus, no arrays larger than 2^32 bytes are possible
 *  1 ) format="appended", all the binary data is at the end <AppendedData encoding="raw"> ... </AppendedData> 
 *  2 ) The first entry after AppendedData tag must be the ascii symbol for underscore
 *  3 ) Length of DataArray in bytes --> sizeof() 
 *  4 ) Actual binary data 
 *  
 *  --> this is adaptable towards parallel output as pvtu (but not HDF5 compatible)
 *  
 *  replacing: plot_cell_data_vtk_writer() with 
 *    - plot_cell_meta_data_vtk_writer
 *    - plot_cell_binary_data_vtk_writer
 *  */
