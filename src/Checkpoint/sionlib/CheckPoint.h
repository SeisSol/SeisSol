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
//#include "Initializer/typedefs.hpp"
#include <iostream>
#include <cstring>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

#include "utils/logger.h"
#include "utils/stringutils.h"
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
	
	/** Identifiers of the files */


      public:
	CheckPoint(unsigned long identifier)
	  : m_identifier(identifier)
	{
	}	

	virtual ~CheckPoint() {}
	
	void setFilename(const char* filename)
	{//initFilename(filename, 0L);
	  initFilename(filename, "sion");
	}
	
	void initLate()
	{
	  int globalrank,numFiles;
	  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE", 0);
	  char fname[1023], *newfname=NULL;
	  m_gComm=comm(); m_lComm = m_gComm; globalrank = rank(); numFiles = 0;
	  seissol::checkpoint::CheckPoint::initLate();
	  	  
#ifdef USE_SIONLIB
	  for (unsigned int i = 0; i < 2; i++) {
	    m_files[i] = sion_paropen_mpi(const_cast<char*>(dataFile(i).c_str()), "bw", &numFiles, m_gComm, &m_lComm,
					  &m_chunksize_sion, &fsblksize, &globalrank, &m_fptr_sion[i], &newfname);
	    fgetpos(m_fptr_sion[i],&m_chunk_pos);
	    checkErr(m_files[i]);
	    checkErr(sion_parclose_mpi(m_files[i]));
	  }
#endif
	}

	void close()
	{

	  //checkErr(::close(m_files[i]));
	  //for (unsigned int i = 0; i < 2; i++){
	    //checkErr(::close(m_files[i]));
	    //logInfo(rank())<<"sionlib:close:sion_parclose_mpi";
	    //logWarning(rank())<<"sionlib:close:sion_parclose_mpi";
	    //checkErr(sion_parclose_mpi(m_files[i]));
	    //logInfo(rank())<<"sionlib:close:sion_parclose_mpi ... done";
	    //logWarning(rank())<<"sionlib:close:sion_parclose_mpi ... done";
	  //}
	}
	
      protected:
	
	int         m_files[2];
	MPI_Comm    m_gComm, m_lComm;
	FILE       *m_fptr_sion[2];
	sion_int64  m_chunksize_sion;
	FILE       *m_fptr;
	fpos_t      m_chunk_pos;
	
	bool exists()
	{
	  if (!seissol::checkpoint::CheckPoint::exists())
	    return false;
	  int file = open();
	  if (file < 0)
	    return false;	  
	  int hasCheckpoint = validate(file);
	  ::close(file);
	  return hasCheckpoint;
	}
	
	/**
	 * Finalize checkpoint writing: Flush the file
	 */
	void finalizeCheckpoint()
	{
	  checkErr(fsync(m_files[odd()]));
	}
	
	/**
	 * Open a check point file
	 *
	 * @return The file handle
	 */
	int open()
	{
	  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE", 0);
	  int globalrank,numFiles;
	  FILE* fptr_sion;
	  char fname[1023], *newfname=NULL;
	  int fh;
	  fh=-1;
	  globalrank = rank(); numFiles = 0; m_gComm = comm(); m_lComm = m_gComm;
#ifdef USE_SIONLIB
	  fh = sion_paropen_mpi(const_cast<char*>(linkFile().c_str()), "br", &numFiles, m_gComm, &m_lComm,
	  			&m_chunksize_sion, &fsblksize, &globalrank, &m_fptr, &newfname);
#endif
	  if (fh < 0)
	    logWarning() << "checkpoint::sionlib::Open():Could not open checkpoint file";
	  return fh;
	}
		
	string linkFile() const
	{
	  string file = std::string(seissol::checkpoint::CheckPoint::linkFile());
	  return file;
	}
		
	string dataFile(int odd) const
	{
	  return seissol::checkpoint::CheckPoint::dataFile(odd); 
	}
	
	/**
	 * @return The current file handle
	 */
	int file() const
	{
	  return m_files[odd()];
	}
	FILE* fptr_sion()
	{
	  return m_fptr_sion[odd()];
	}
	FILE* fptr()
	{
	  return m_fptr;
	}
	fpos_t chunk_pos() const
	{
	  return m_chunk_pos;
	}
	void set_chunksize(sion_int64 chunksize)
	{
	  m_chunksize_sion = chunksize;
	}
	sion_int64 chunksize()
	{
	  return m_chunksize_sion;
	}	
	/**
	 * @return The identifier of the file
	 */
	unsigned long identifier() const
	{
	  return m_identifier;
	}

      private:
	/**
	 * Validate an existing check point file
	 */
	bool validate(int file) const
	{
	  unsigned long id;
	  checkErr(sion_fread(&id, sizeof(unsigned long),1,file));
	  if (id != identifier()) {
	    logWarning() << "Checkpoint identifier does not match";
	    return false;
	  }
	  return true;
	}
	
      protected:
	template<typename T>
	static void checkErr(T ret)
	{
	  if (ret < 0)
	    logError() << "Error in the SIONlib checkpoint module:" << strerror(errno);
	}
	
	/**
	 * Can be used to check read/write errors
	 *
	 * @param ret The return value
	 * @param target The expected return value (> 0)
	 */
	template<typename T, typename U>
	static void checkErr(T ret, U target)
	{
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
