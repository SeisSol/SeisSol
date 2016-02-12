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
 * Main class for SIONlib checkpoints
 */

#ifndef CHECKPOINT_SIONLIB_CHECK_POINT_H
#define CHECKPOINT_SIONLIB_CHECK_POINT_H

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#ifdef USE_SIONLIB
#include <sion.h>
#else
#define sion_int32 int
#define sion_int64 long long
#endif
#include <string>
#include <cstring>
#include <errno.h>
#include <sys/stat.h>
#include "utils/logger.h"
#include "utils/env.h"
#include "Checkpoint/CheckPoint.h"
#include "Initializer/preProcessorMacros.fpp"

using namespace std;

namespace seissol
{  
  namespace checkpoint
  {   
    namespace sionlib
    {   
      class CheckPoint : virtual public seissol::checkpoint::CheckPoint
      {
      private:
	/** Checkpoint identifier (written to the beginning of the file) */
	const unsigned long m_identifier;
	
      public:
	CheckPoint(unsigned long identifier)
	  : m_identifier(identifier){}	

	virtual ~CheckPoint() {}
	
	void setFilename(const char* filename) {
	  initFilename(filename, "sion");
	}
	
	void initLate() {
#ifdef USE_SIONLIB
 	  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE", 0);
	  int globalrank,numFiles;
	  char fname[1023], *newfname=NULL;
	  m_gComm=comm(); m_lComm = m_gComm; globalrank = rank(); numFiles = 0;
	  seissol::checkpoint::CheckPoint::initLate();
	  	  
	  for (unsigned int i = 0; i < 2; i++) {
	    m_files[i] = sion_paropen_mpi(const_cast<char*>(dataFile(i).c_str()), "bw", &numFiles, m_gComm, &m_lComm,
					  &m_chunksize, &fsblksize, &globalrank, &m_fptr[i], &newfname);
	    fgetpos(m_fptr[i],&m_chunkpos);
	    checkErr(m_files[i]);
	    logInfo(rank())<<"opened sionlib checkpoint:chunksize"<<m_chunksize<<"\n";
	  }
#endif
	}
	/** close single checkpoint file */
	void close_file(int fh){
	  checkErr(sion_parclose_mpi(fh));
	}
	/** close checkpoint */
	void close(){
	  for (unsigned int i = 0; i < 2; i++){
	    checkErr(sion_parclose_mpi(m_files[i]));
	  }
	}
	
      protected:
	int         m_files[2];
	MPI_Comm    m_gComm, m_lComm;
	FILE       *m_fptr[2];
	sion_int64  m_chunksize;
	fpos_t      m_chunkpos;
	
	bool exists() {
	  if (!seissol::checkpoint::CheckPoint::exists())
	    return false;
	  int file = open();
	  if (file < 0)
	    return false;	  
	  int hasCheckpoint = validate(file);
	  close_file(file);//::close(file);
	  return hasCheckpoint;
	}
		
	void flushCheckpoint() {	
		SCOREP_USER_REGION_DEFINE(r_flush);
		SCOREP_USER_REGION_BEGIN(r_flush, "checkpoint_flush", SCOREP_USER_REGION_TYPE_COMMON);
		checkErr(fflush(m_fptr[odd()]));
		SCOREP_USER_REGION_END(r_flush);
	}
	
	/** Open a check point file, @return The file handle */
	int open() {
	  int fh; fh=-1;
#ifdef USE_SIONLIB
	  int globalrank,numFiles; FILE* fptr;
	  char fname[1023], *newfname=NULL;
	  globalrank = rank(); numFiles = 0; m_gComm = comm(); m_lComm = m_gComm;
	  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE", 0);
	  fh = sion_paropen_mpi(const_cast<char*>(linkFile().c_str()), "br", &numFiles, m_gComm, &m_lComm,
	  			&m_chunksize, &fsblksize, &globalrank, &fptr, &newfname);
#endif
	  if (fh < 0)
	    logWarning() << "checkpoint::sionlib::Open():Could not open checkpoint file";
	  return fh;
	}
		
	string linkFile() const {
	  string file = std::string(seissol::checkpoint::CheckPoint::linkFile());
	  return file;
	}
		
	string dataFile(int odd) const {
	  return seissol::checkpoint::CheckPoint::dataFile(odd); 
	}	
	/** @return The current file handle */
	int file() const {
	  return m_files[odd()];
	}
	/** @return The current file pointer */
	FILE* fptr() {
	  return m_fptr[odd()];
	}
	/** reset file position to beginn of the chunk */
	void setpos(){
	  fsetpos(m_fptr[odd()],&m_chunkpos);
	}
	/** set chunk size */
	void set_chunksize(sion_int64 chunksize) {
	  m_chunksize = chunksize;
	}
	/** @return The identifier of the file */
	unsigned long identifier() const {
	  return m_identifier;
	}

      private:

	/** Validate an existing check point file */
	bool validate(int file) const
	{
	  unsigned long id;
#ifdef USE_SIONLIB
	  checkErr(sion_fread(&id, sizeof(unsigned long),1,file));
#endif
	  if (id != identifier()) {
	    logWarning() << "Checkpoint identifier does not match";
	    return false;
	  }
	  return true;
	}
	
      protected:
	template<typename T>
	static void checkErr(T ret) {
	  if (ret < 0)
	    logError() << "Error in the SIONlib checkpoint module:" << strerror(errno);
	}	
	/** Can be used to check read/write errors
	 * @param ret The return value
	 * @param target The expected return value (> 0)
	 */
	template<typename T, typename U>
	static void checkErr(T ret, U target) {
	  checkErr(ret);
	  if (ret != target)
	    logError() << "Error in the SIONlib checkpoint module:"
		       << target << "bytes expected;" << ret << "bytes gotten";
	}
      }; 
    }
  }
}
#endif // CHECKPOINT_SIONLIB_CHECK_POINT_H
