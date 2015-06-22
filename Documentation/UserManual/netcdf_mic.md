# Building netCDF for native Xeon Phi execution
## Information
 * contact: Alex Breuer (breuer@mytum.de)
 * This guide was tested on the MAC-Cluster of TUM (14/12/22)
 * Yes it's a mess..
 * No efforts have been made to generalize the instructions; please adopt to your settings!

## zlib
* get zlib
  wget http://downloads.sourceforge.net/project/libpng/zlib/1.2.8/zlib-1.2.8.tar.gz

* configure zlib
  ./configure --prefix=/home/hpc/pr63so/di56dok/software/netcdf/mic/zlib-1.2.8/build

* edit Makefile and add "-mmic"s and MIC-libraries
  CFLAGS=-mmic -O3  -D_LARGEFILE64_SOURCE=1 -DHAVE_HIDDEN
  LDFLAGS=-mmic
  SFLAGS=-mmic -O3  -fPIC -D_LARGEFILE64_SOURCE=1 -DHAVE_HIDDEN -mmic
  examplesh$(EXE): example.o $(SHAREDLIBV)
          $(CC) $(CFLAGS) -o $@ example.o -L. $(SHAREDLIBV)  -L/lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/compiler/lib/mic/ -limf -lsvml -lirng -lintlc

  minigzipsh$(EXE): minigzip.o $(SHAREDLIBV)
          $(CC) $(CFLAGS) -o $@ minigzip.o -L. $(SHAREDLIBV)  -L/lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/compiler/lib/mic/ -limf -lsvml -lirng -lintlc

  example64$(EXE): example64.o $(STATICLIB)
          $(CC) $(CFLAGS) -o $@ example64.o $(TEST_LDFLAGS) -L/lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/compiler/lib/mic/ -limf -lsvml -lirng -lintlc

  minigzip64$(EXE): minigzip64.o $(STATICLIB)
          $(CC) $(CFLAGS) -o $@ minigzip64.o $(TEST_LDFLAGS) -L/lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/compiler/lib/mic/ -limf -lsvml -lirng -lintlc

* make
  make
  make install

## HDF5
* get hdf5
  wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.14.tar

* configure hdf5
  ./configure CC=mpiicc CXX=mpiicpc FC=mpiifort --with-zlib=/home/hpc/pr63so/di56dok/software/netcdf/mic/zlib-1.2.8/build  --host=x86_64-unknown-linux --prefix=/home/hpc/pr63so/di56dok/software/netcdf/mic/hdf5-1.8.14/build --enable-parallel CFLAGS='-mt_mpi' LD_FLAGS='-mt_mpi' FCFLAGS='-mt_mpi' CXXFLAGS='-mt_mpi' --enable-shared=no

* duplicate the HDF5-folder and build a host only version
  cp -r hdf5-1.8.14 hdf5-1.8.14_host
  make -j

* (at the moment the upcoming build fails) copy over the required binaries from the host-version to the mic-version and continue building
  cp ../hdf5-1.8.14_host/src/H5make_libsettings src/
  cp ../hdf5-1.8.14_host/src/H5detect src/
  

* add "-mmic"s and MIC libraries to config.status
  S["CXXFLAGS"]="-mmic -mt_mpi"
  S["FCFLAGS_f90"]="-mmic -mt_mpi"
  S["FCFLAGS"]="-mmic -mt_mpi"
  S["LDFLAGS"]="-mmic -mt_mpi -L/lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/lib/mic/ -limf -lsvml -lirng -lintlc"
  S["CFLAGS"]="-mmic -mt_mpi"
  S["AM_LDFLAGS"]=" -mmic -mt_mpi -L/home/hpc/pr63so/di56dok/software/netcdf/mic/hdf5-1.8.14/hdf5/lib -L/home/hpc/pr63so/di56dok/software/netcdf/mic/zlib-1.2.8/build/lib -L/lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/compiler/lib/mic/ -limf -lsvml -lirng -lintlc"
  S["AM_CPPFLAGS"]=" -mmic -mt_mpi -I/home/hpc/pr63so/di56dok/software/netcdf/mic/zlib-1.2.8/build/include -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE "
  S["AM_CXXFLAGS"]="-mmic -mt_mpi"
  S["AM_FCFLAGS"]="-mmic -mt_mpi"
  S["AM_CFLAGS"]="-mmic -mt_mpi"
  S["H5_LDFLAGS"]="-mmic -mt_mpi -L/lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/compiler/lib/mic/ -limf -lsvml -lirng -lintlc"
  S["H5_CXXFLAGS"]=" -mmic -mt_mpi"
  S["H5_FCFLAGS"]=" -mmic -mt_mpi -O3"
  S["H5_CPPFLAGS"]=" -mmic -mt_mpi -D_POSIX_C_SOURCE=199506L   -DNDEBUG -UH5_DEBUG_API"
  S["H5_CFLAGS"]="-mmic -mt_mpi -std=c99  -O3"
  S["CPPFLAGS"]="-mmic -mt_mpi -DpgiFortran"

* build the library
  make
  make install

## netCDF
* get netcdf
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.2.tar.gz

* configure netcdf
./configure CC=mpiicc CXX=mpiicpc FC=mpiifort CFLAGS="-mt_mpi -mmic -I/home/hpc/pr63so/di56dok/software/netcdf/mic/hdf5-1.8.14/build/include" LDFLAGS="-mmic -mt_mpi -L/home/hpc/pr63so/di56dok/software/netcdf/mic/hdf5-1.8.14/build/lib -L/home/hpc/pr63so/di56dok/software/netcdf/mic/zlib-1.2.8/build/lib -L/lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/compiler/lib/mic/" LIBS="-limf -lsvml -lirng -lintlc" --host=x86_64-unknown-linux --prefix=/home/hpc/pr63so/di56dok/software/netcdf/mic/netcdf-4.3.2/build --enable-shared=no

* patch netcdf according to https://github.com/Unidata/netcdf-c/blob/435d8a03ed28bb5ad63aff12cbc6ab91531b6bc8/libsrc4/nc4file.c
  cd libsrc4/
  mv nc4file.c nc4file.c_back
  wget https://github.com/Unidata/netcdf-c/raw/435d8a03ed28bb5ad63aff12cbc6ab91531b6bc8/libsrc4/nc4file.c --no-check-certificate
* build netcdf
  make -j
  make install

* OVER finally!
  +-------------------------------------------------------------+
  | Congratulations! You have successfully installed netCDF!    |
  [...]

## SeisSol
* set the libraries
  export LIBRARY_PATH=/lrz/sys/intel/compiler140_106/composer_xe_2013_sp1.1.106/compiler/lib/mic/:$LIBRARY_PATH
  export LIBRARY_PATH=/home/hpc/pr63so/di56dok/software/netcdf/mic/zlib-1.2.8/build/lib/:$LIBRARY_PATH
  export LD_LIBRARY_PATH=/home/hpc/pr63so/di56dok/software/netcdf/mic/zlib-1.2.8/build/lib/:$LD_LIBRARY_PATH
* and in Scons
  netcdf=yes hdf5Dir=/home/hpc/pr63so/di56dok/software/netcdf/mic/hdf5-1.8.14/build/ netcdfDir=/home/hpc/pr63so/di56dok/software/netcdf/mic/netcdf-4.3.2/build/
