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
//#include "sionIOcom.h"
#include <string>
#include <cstring>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include "utils/logger.h"
#include "utils/env.h"
#include "Checkpoint/CheckPoint.h"
#include "Initializer/preProcessorMacros.fpp"

using namespace std;

class SionIOGroup {
  
private:

protected:
  
#ifdef USE_MPI
  MPI_Comm newcomm;
#endif //USE_MPI
  int group;
  int key;
  int numfiles;
  
public: 
  
  SionIOGroup(void) : group(1),key(0),numfiles(1) {}
  SionIOGroup(MPI_Comm aComm){
    string val,substr,subn;
    int rank,size;
    
    try{
      val = string(getenv("SEISSOL_CHECKPOINT_SION_NUM_FILES"));
      numfiles = atoi(val.c_str());
      rank=0;size=1;
#ifdef USE_MPI
      MPI_Comm_rank(aComm,&rank);
      MPI_Comm_size(aComm,&size);
#endif //USE_MPI
      group = (rank/(size/numfiles))+1;
    } 
    catch(...){
      try{
	val = string(getenv("SEISSOL_CHECKPOINT_SION_GROUP"));
	if (val==string("HOSTNAME")){
	  substr = string(getenv("SEISSOL_CHECKPOINT_SION_SUBSTR"));
	  subn = string(getenv("SEISSOL_CHECKPOINT_SION_SUBN"));
	  set_group(substr,subn);
	  cout<<"SEISSOL_CHECKPOINT_SION_GROUP:HOSTNAME"<<endl;
	}else{group = atoi(val.c_str());}}
      catch(...){group=1;cout<<"SEISSOL_CHECKPOINT_SION_GROUP:fallback"<<endl;}}
#ifdef USE_MPI
    MPI_Comm_split(aComm,group,key,&newcomm);
    logInfo(rank)<<"|group:"<<group<<"|key:"<<key<<"|numfiles:"<<numfiles;
#endif //USE_MPI
  }
  int get_group()  {return group  ;}
#ifdef USE_MPI
  int get_newcomm(){return newcomm;}
#endif //USE_MPI
  int get_key()    {return key    ;}  

private:

  
  void set_group(string substr,string subn){int n;
    try{  n = atoi(subn.c_str());
      group = get_from_substr(hostname(),substr,n);}
    catch(...){group = 1;}
    return;
  }
  string hostname(){ char name2[256]; int iret;
    iret = gethostname (name2, sizeof name2);
    return string(name2);
  }
  int get_from_substr(string name,string sub, int n=2){
    int iret;  size_t keyloc; string dum;
    keyloc = name.find(sub);
    dum = name.substr(keyloc+1,keyloc+n);
    return atoi(dum.c_str());
  }
};

namespace seissol{  

  namespace checkpoint{   

    namespace sionlib{
      
      class CheckPoint : virtual public seissol::checkpoint::CheckPoint{

      private:
	/** Checkpoint identifier (written to the beginning of the file) */
	const unsigned long m_identifier;

      public:

	CheckPoint(unsigned long identifier)
	  : m_identifier(identifier),m_flush(utils::Env::get<int>("SEISSOL_CHECKPOINT_SION_FLUSH",0)){}	
	
	virtual ~CheckPoint() {}
	//~CheckPoint();
	
	void setFilename(const char* filename) {initFilename(filename, 0L);}
	
	void initLate() {

#ifdef USE_SIONLIB
	  seissol::checkpoint::CheckPoint::initLate();
	  // Create the folder
	  if (rank() == 0) {
	    for (int i = 0; i < 2; i++) {
	      int ret = mkdir(seissol::checkpoint::CheckPoint::dataFile(i).c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
	      if (ret < 0 && errno != EEXIST){ checkErr(ret); }
	    }
	  }	  
#ifdef USE_MPI // Make sure all processes see the folders
	  MPI_Barrier(comm());
#endif // USE_MPI	  
	  writeinit();
#endif
	}
	virtual void writeinit(){}
	
	/** close single checkpoint file */
	void close_file(int fh){
	  checkErr(sion_parclose_mpi(fh));
	  logInfo(rank())<<"close_file()";
	}
	/** close checkpoint */
	void close(){
	  for (unsigned int i = 0; i < 2; i++){
	    checkErr(sion_parclose_mpi(m_files[i]));
	  }
	  logInfo(rank())<<"close()";
	}
	
      protected:
		    
	int         m_files[2];
#ifdef USE_MPI // Make sure all processes see the folders
	MPI_Comm    m_gComm, m_lComm;
#endif //USE_MPI // Make sure all processes see the folders
	FILE       *m_fptr[2];
	sion_int64  m_chunksize;
	fpos_t      m_chunkpos;
	SionIOGroup m_iogroup;
	int         m_flush;
		    
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
	    if(m_flush==0){
	      checkErr(sion_flush(m_files[odd()]));	      
	    }else{
	      checkErr(fflush(m_fptr[odd()]));
	    }
	    SCOREP_USER_REGION_END(r_flush);
	    logInfo(rank())<<"flush_Checkpoint()";
	  }
	
	void open_all() {	  
	  for (unsigned int i = 0; i < 2; i++) {
	    m_files[i] = open();
	  }
	  return ;
	}
	/** Open a check point file, @return The file handle */
	int open() {
	  int fh; fh=-1;
#ifdef USE_SIONLIB
	  int globalrank,numFiles; FILE* fptr;
	  char fname[1023], *newfname=NULL;
	  globalrank = rank(); numFiles = -1; m_gComm = comm(); m_lComm = m_iogroup.get_newcomm();
	  sion_int32 fsblksize= utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE",-1);
	  logInfo(rank())<<":: open:file:"<<linkFile().c_str()<<"|group:"<<m_iogroup.get_group();
	  fh = sion_paropen_mpi(const_cast<char*>(linkFile().c_str()), "br", &numFiles, m_gComm, &m_lComm,
	  			&m_chunksize, &fsblksize, &globalrank, &fptr, &newfname);
#endif
	  if (fh < 0){logWarning() << "checkpoint::sionlib::Open():Could not open checkpoint file";}
	  return fh;
	}
		
	string linkFile(){return std::string(seissol::checkpoint::CheckPoint::linkFile())+ "/checkpoint";}
		
	string dataFile(int odd) {return seissol::checkpoint::CheckPoint::dataFile(odd) + "/checkpoint";}	
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
