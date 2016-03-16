//-*-c++-*-
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Gilbert Brietzke (gilbert.brietzke AT lrz.de, http://www.lrz.de)
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
#include "Fault.h"
#include "Kernels/precision.hpp"

bool seissol::checkpoint::sionlib::Fault::init(double* mu, double* slipRate1, double* slipRate2,
					       double* slip, double* slip1,double* slip2,
					       double* state, double* strength,
					       unsigned int numSides, unsigned int numBndGP){
#ifdef USE_SIONLIB
  seissol::checkpoint::Fault::init(mu, slipRate1, slipRate2, slip, slip1, slip2, state, strength, numSides, numBndGP);
  set_chunksize(this->numSides()*this->numBndGP()*sizeof(real) + sizeof(int) + sizeof(unsigned long));
  if (numSides == 0)
    return true;
  m_iogroup = SionIOGroup(comm());
#endif
  return exists();
}

void seissol::checkpoint::sionlib::Fault::writeinit(){
#ifdef USE_SIONLIB
  int globalrank,numFiles;
  char *newfname=NULL;
  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE",-1);  
  unsigned long lidentifier;
  int timestepFault;
  
  m_gComm = comm(); m_lComm = m_iogroup.get_newcomm();
  globalrank = rank(); numFiles = -1;
  lidentifier = identifier();timestepFault=0;
  for (unsigned int i = 0; i < 2; i++) {
    logInfo(rank())<<"writeinit: connect to file:"<<dataFile(i).c_str()<<"|group:"<<m_iogroup.get_group();
    m_files[i] = sion_paropen_mpi(const_cast<char*>(dataFile(i).c_str()), "bw", &numFiles, m_gComm, &m_lComm,
				  &m_chunksize, &fsblksize, &globalrank, &m_fptr[i], &newfname);
    fgetpos(m_fptr[i],&m_chunkpos);//get position of datachunk for subsequent writes
    checkErr(m_files[i]);
    checkErr(sion_fwrite(&lidentifier, sizeof(unsigned long),1,m_files[i]));
    checkErr(sion_fwrite(&timestepFault, sizeof(timestepFault),1,m_files[i]));    
    if (numSides() > 0){
      for (int j = 0; j < NUM_VARIABLES; j++){
	checkErr(sion_fwrite(data(j),sizeof(real),this->numSides()*this->numBndGP(),m_files[i]));
      }
    }
    close_file(m_files[i]);
  }
  for (unsigned int i = 0; i < 2; i++) {
    logInfo(rank())<<"writeinit: reconnect to file:"<<dataFile(i).c_str()<<"|group:"<<m_iogroup.get_group();
    m_files[i] = sion_paropen_mpi(const_cast<char*>(dataFile(i).c_str()), "bw", &numFiles, m_gComm, &m_lComm,
				  &m_chunksize, &fsblksize, &globalrank, &m_fptr[i], &newfname);
  }
#endif
}

void seissol::checkpoint::sionlib::Fault::load(int &timestepFault) {
#ifdef USE_SIONLIB
  if (numSides() == 0)
    return;  
  logInfo(rank()) << "Loading fault checkpoint";
  seissol::checkpoint::CheckPoint::load();
  int file; FILE *file_ptr; char *newfname=NULL;
  int globalrank,numFiles;  unsigned long lidentifier;  
  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE",-1);
  m_gComm = comm(); m_lComm = m_iogroup.get_newcomm();
  globalrank  = rank(); numFiles = -1;
  file = sion_paropen_mpi(const_cast<char*>(linkFile().c_str()), "br", &numFiles, m_gComm, &m_lComm,
			  &m_chunksize, &fsblksize, &globalrank, &file_ptr, &newfname);
  checkErr(file);
  checkErr(sion_fread(&lidentifier, sizeof(unsigned long),1,file));
  checkErr(sion_fread(&timestepFault, sizeof(timestepFault),1,file));
  if (numSides() > 0){
    for (int i = 0; i < NUM_VARIABLES; i++)
      checkErr(sion_fread(data(i),sizeof(real),this->numSides()*this->numBndGP(),file));
  }
  if (ferror (file_ptr))
    logWarning(rank())<<"Error reading fault data SIONlib-checkpoint\n";
  checkErr(sion_parclose_mpi(file));
  logInfo(rank()) << "Loading fault SIONlib-checkpoint done";
#endif
}

void seissol::checkpoint::sionlib::Fault::write(int timestepFault) {  
  SCOREP_USER_REGION("CheckPointFault_write", SCOREP_USER_REGION_TYPE_FUNCTION);
#ifdef USE_SIONLIB
  if (numSides() == 0)
    return;
  logInfo(rank()) << "Writing fault check point.";
  if (m_method == 0){
    open_write_close(timestepFault);
  }else{
    
    unsigned long lidentifier;
    lidentifier = identifier();
    setpos();
    
    SCOREP_USER_REGION_DEFINE(r_write_header);
    SCOREP_USER_REGION_BEGIN(r_write_header, "checkpoint_write_fault_header", SCOREP_USER_REGION_TYPE_COMMON);
    
    checkErr(sion_fwrite(&lidentifier, sizeof(unsigned long),1,m_files[odd()]));
    checkErr(sion_fwrite(&timestepFault, sizeof(timestepFault),1,m_files[odd()]));
    
    SCOREP_USER_REGION_END(r_write_header);
    
    SCOREP_USER_REGION_DEFINE(r_write_fault);
    SCOREP_USER_REGION_BEGIN(r_write_fault, "checkpoint_write_fault", SCOREP_USER_REGION_TYPE_COMMON);
    if (numSides() > 0){
      for (int i = 0; i < NUM_VARIABLES; i++){
	checkErr(sion_fwrite(data(i),sizeof(real),this->numSides()*this->numBndGP(),m_files[odd()]));
      }
    }
    SCOREP_USER_REGION_END(r_write_fault);
    flushCheckpoint(); 
    //  finalizeCheckpoint();  
  }
  logInfo(rank()) << "Writing fault check point. Done.";
#endif
}

void seissol::checkpoint::sionlib::Fault::open_write_close(int timestepFault){
#ifdef USE_SIONLIB
  int globalrank,numFiles;
  char *newfname=NULL;
  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE",-1);  
  unsigned long lidentifier;
  lidentifier = identifier();
  m_gComm = comm(); m_lComm = m_iogroup.get_newcomm();
  globalrank = rank(); numFiles = -1;
  
  //logInfo(rank())<<"opening file for write:"<<dataFile(odd()).c_str();
  m_files[odd()] = sion_paropen_mpi(const_cast<char*>(dataFile(odd()).c_str()), "bw", &numFiles, m_gComm, &m_lComm,
				    &m_chunksize, &fsblksize, &globalrank, &m_fptr[odd()], &newfname);
  checkErr(m_files[odd()]);
  //logInfo(rank())<<"opened file for write:"<<linkFile().c_str();
  
  SCOREP_USER_REGION_DEFINE(r_write_header);
  SCOREP_USER_REGION_BEGIN(r_write_header, "checkpoint_write_fault_header", SCOREP_USER_REGION_TYPE_COMMON);
  checkErr(sion_fwrite(&lidentifier, sizeof(unsigned long),1,m_files[odd()]));
  checkErr(sion_fwrite(&timestepFault, sizeof(timestepFault),1,m_files[odd()]));    
  SCOREP_USER_REGION_END(r_write_header);
  
  SCOREP_USER_REGION_DEFINE(r_write_fault);
  SCOREP_USER_REGION_BEGIN(r_write_fault, "checkpoint_write_fault", SCOREP_USER_REGION_TYPE_COMMON);  
  //  if (numSides() > 0){
  for (int j = 0; j < NUM_VARIABLES; j++){
    checkErr(sion_fwrite(data(odd()),sizeof(real),this->numSides()*this->numBndGP(),m_files[odd()]));
  }
  //}
  SCOREP_USER_REGION_END(r_write_fault);

  checkErr(sion_parclose_mpi(m_files[odd()]));
  //logInfo(rank())<<"closed file for write:"<<dataFile(odd()).c_str();
#endif
}
