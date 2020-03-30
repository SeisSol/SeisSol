# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *

class Pumgen(SConsPackage):
    homepage = "https://github.com/SeisSol/PUMGen/wiki/How-to-compile-PUMGen"
    version('develop',
            git='https://github.com/SeisSol/PUMGen.git',
            branch='master',
            submodules=True)

    maintainers = ['ravil-mobile']
    variant('simmetrix_support', default=False)
    depends_on('mpi')
        
    depends_on('netcdf-c +shared +mpi') # NOTE: only tested with 4.4.0 version
    depends_on('hdf5 +fortran +shared +mpi') # NOTE: only tested with 1.8.21 version
    depends_on('pumi +int64 +zoltan -fortran', when='~simmetrix_support')
    depends_on('pumi +int64 simmodsuite=kernels +zoltan -fortran', when='+simmetrix_support')
    depends_on('zoltan@3.83 +parmetis+int64 -fortran')

    def build_args(self, spec, prefix):                                                                               
        args=[]                                                                                                                                                                                                                         
        mpi_id = spec['mpi'].name + spec['mpi'].version.string                                                        
        args.append('mpiLib=' + mpi_id)                                                                               
        args.append('cc=mpicc')                                                                                       
        args.append('cxx=mpicxx')                                                                                                                                                                                                   
        if '+simmetrix_support' in spec:     
            args.append('simModSuite=yes')                                                                            
        return args                                                                                                   
                                                                                                                  
    def build(self, spec, prefix):                                                                                    
        args = self.build_args(spec, prefix)                                                                          
        scons(*args)                                                                                                  
    
    def install(self,spec,prefix):
        install_tree('build',prefix.bin)
