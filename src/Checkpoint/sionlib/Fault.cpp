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

bool seissol::checkpoint::sionlib::Fault::init(double* mu, double* slipRate1, double* slipRate2, double* slip, double* slip1,
					       double* slip2, double* state, double* strength,
					       unsigned int numSides, unsigned int numBndGP){
#ifdef USE_SIONLIB
  seissol::checkpoint::Fault::init(mu, slipRate1, slipRate2, slip, slip1, slip2, state, strength, numSides, numBndGP);
  if (numSides == 0)
    return true;
  set_chunksize(this->numSides()*this->numBndGP()*sizeof(real) + sizeof(int) + sizeof(unsigned long));
#endif
  return exists();
}

void seissol::checkpoint::sionlib::Fault::load(int &timestepFault) {
#ifdef USE_SIONLIB
  if (numSides() == 0)
    return;  
  seissol::checkpoint::CheckPoint::load();
  int file; FILE *file_ptr; char fname[1023], *newfname=NULL;
  int globalrank,numFiles;  unsigned long lidentifier;  
  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE", 0);
  m_gComm = comm(); m_lComm = m_gComm;	
  globalrank  = rank(); numFiles = 0;
  file = sion_paropen_mpi(const_cast<char*>(linkFile().c_str()), "br", &numFiles, m_gComm, &m_lComm,
			  &m_chunksize, &fsblksize, &globalrank, &file_ptr, &newfname);
  checkErr(file);
  checkErr(sion_fread(&lidentifier, sizeof(unsigned long),1,file));
  checkErr(sion_fread(&timestepFault, sizeof(timestepFault),1,file));
  for (int i = 0; i < NUM_VARIABLES; i++)
    checkErr(sion_fread(data(i),sizeof(real),this->numSides()*this->numBndGP(),file));
  if (ferror (file_ptr))
    logWarning(rank())<<"Error reading fault data SIONlib-checkpoint\n";
  checkErr(sion_parclose_mpi(file));
#endif
}

void seissol::checkpoint::sionlib::Fault::write(int timestepFault) {  
#ifdef USE_SIONLIB
  if (numSides() == 0)
    return;
  int globalrank,numFiles;
  char fname[1023], *newfname=NULL; 
  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE", 0);
  unsigned long lidentifier;
  lidentifier = identifier();
  m_gComm = comm(); m_lComm = m_gComm;
  globalrank = rank(); numFiles = 0; 
  setpos();
  checkErr(sion_fwrite(&lidentifier, sizeof(unsigned long),1,m_files[odd()]));
  checkErr(sion_fwrite(&timestepFault, sizeof(timestepFault),1,m_files[odd()]));
  for (int i = 0; i < NUM_VARIABLES; i++){
    checkErr(sion_fwrite(data(i),sizeof(real),this->numSides()*this->numBndGP(),m_files[odd()]));
  }
  flushCheckpoint(); 
  //  finalizeCheckpoint();  
#endif
}
