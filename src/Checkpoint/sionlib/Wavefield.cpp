//*-* c++ -*-
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Gilbert Brietzke (gilbert.brietzke AT LRZ.de, http://www.lrz.de)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 *
 * @section DESCRIPTION
 */
#include "Wavefield.h"


bool seissol::checkpoint::sionlib::Wavefield::init(real* dofs, unsigned int numDofs){
  seissol::checkpoint::Wavefield::init(dofs, numDofs);
#ifdef USE_SIONLIB
  set_chunksize(this->numDofs()*sizeof(double)+sizeof(double)+sizeof(int)+sizeof(unsigned long));
#endif
  return exists();
}

void seissol::checkpoint::sionlib::Wavefield::writeinit(){
#ifdef USE_SIONLIB
  int globalrank,numFiles;
  char fname[1023], *newfname=NULL;
  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE",-1);  
  unsigned long lidentifier;
  double time;int timestepWaveField;
  
  m_gComm = comm(); m_lComm = m_iogroup.get_newcomm();
  globalrank = rank(); numFiles = -1;
  
  for (unsigned int i = 0; i < 2; i++) {
    logInfo(rank())<<"writeinit: connect to file:"<<dataFile(i).c_str()<<"|group:"<<m_iogroup.get_group()<<"|lcomm:"<<m_lComm;
    m_files[i] = sion_paropen_mpi(const_cast<char*>(dataFile(i).c_str()), "bw", &numFiles, m_gComm, &m_lComm,
				  &m_chunksize, &fsblksize, &globalrank, &m_fptr[i], &newfname);
    fgetpos(m_fptr[i],&m_chunkpos);
    checkErr(m_files[i]);
    checkErr(sion_fwrite(&lidentifier, sizeof(unsigned long),1,m_files[odd()]));
    checkErr(sion_fwrite(&time, sizeof(time),1,m_files[odd()]));
    checkErr(sion_fwrite(&timestepWaveField, sizeof(timestepWaveField),1,m_files[odd()]));
    checkErr(sion_fwrite(dofs(),sizeof(real), numDofs(),m_files[odd()]));  
    ::close(m_files[i]);
    m_files[i] = sion_paropen_mpi(const_cast<char*>(dataFile(i).c_str()), "bw", &numFiles, m_gComm, &m_lComm,
				  &m_chunksize, &fsblksize, &globalrank, &m_fptr[i], &newfname);
    checkErr(m_files[i]);
  }
  for (unsigned int i = 0; i < 2; i++) {
    logInfo(rank())<<"writeinit: reconnect to file:"<<dataFile(i).c_str()<<"|group:"<<m_iogroup.get_group()<<"|lcomm:"<<m_lComm;
    m_files[i] = sion_paropen_mpi(const_cast<char*>(dataFile(i).c_str()), "bw", &numFiles, m_gComm, &m_lComm,
				  &m_chunksize, &fsblksize, &globalrank, &m_fptr[i], &newfname);
  }
#endif
}

void seissol::checkpoint::sionlib::Wavefield::load(double &time, int &timestepWaveField) {  
#ifdef USE_SIONLIB
  logInfo(rank()) << "Loading wave field checkpoint";
  seissol::checkpoint::CheckPoint::load();
  int globalrank,numFiles,file; FILE *file_ptr;
  char fname[1023], *newfname=NULL;
  unsigned long lidentifier;
  real *dofs1;    
  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE",-1);
  m_gComm = comm(); m_lComm = m_iogroup.get_newcomm();
  globalrank = rank(); numFiles = -1; 
  file = sion_paropen_mpi(const_cast<char*>(linkFile().c_str()), "br", &numFiles, m_gComm, &m_lComm,
			  &m_chunksize, &fsblksize, &globalrank, &file_ptr, &newfname);
  checkErr(file);
  checkErr(sion_fread(&lidentifier, sizeof(unsigned long),1,file));
  checkErr(sion_fread(&time, sizeof(time),1,file));
  checkErr(sion_fread(&timestepWaveField, sizeof(timestepWaveField),1,file));
  checkErr(sion_fread(dofs(), sizeof(double),numDofs(),file));
  if (ferror (file_ptr))
    logWarning(rank())<<"Error reading dofs\n";
  checkErr(sion_parclose_mpi(file));
  logInfo(rank()) << "Loading wave field checkpoint done";
#endif
}

void seissol::checkpoint::sionlib::Wavefield::write(double time, int timestepWaveField) {
#ifdef USE_SIONLIB

  SCOREP_USER_REGION("CheckPoint_write", SCOREP_USER_REGION_TYPE_FUNCTION);
  logInfo(rank()) << "Writing check point.";

  int globalrank,numFiles;
  char fname[1023], *newfname=NULL;
  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE",-1);  
  unsigned long lidentifier;
  //m_gComm = comm(); m_lComm = m_gComm;
  globalrank = rank(); numFiles = -1;
  lidentifier = identifier();
  setpos();

  SCOREP_USER_REGION_DEFINE(r_write_header);
  SCOREP_USER_REGION_BEGIN(r_write_header, "checkpoint_write_header", SCOREP_USER_REGION_TYPE_COMMON);

  checkErr(sion_fwrite(&lidentifier, sizeof(unsigned long),1,m_files[odd()]));
  checkErr(sion_fwrite(&time, sizeof(time),1,m_files[odd()]));
  checkErr(sion_fwrite(&timestepWaveField, sizeof(timestepWaveField),1,m_files[odd()]));

  SCOREP_USER_REGION_END(r_write_header);
  SCOREP_USER_REGION_DEFINE(r_write_wavefield);
  SCOREP_USER_REGION_BEGIN(r_write_wavefield, "checkpoint_write_wavefield", SCOREP_USER_REGION_TYPE_COMMON);

  checkErr(sion_fwrite(dofs(),sizeof(real), numDofs(),m_files[odd()]));  

  SCOREP_USER_REGION_END(r_write_wavefield);

  flushCheckpoint();
  //  finalizeCheckpoint();  
  logInfo(rank()) << "Writing check point. Done.";
#endif
}
