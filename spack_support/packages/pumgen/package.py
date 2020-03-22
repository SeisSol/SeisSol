# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install pumgen
#
# You can edit this file again by typing:
#
#     spack edit pumgen
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *


class Pumgen(SConsPackage):
    homepage = "https://github.com/SeisSol/PUMGen/wiki/How-to-compile-PUMGen"
    version('develop',
            git='https://github.com/SeisSol/PUMGen.git',
            branch='master')

    maintainers = ['ravil-mobile']
    variant('simmetrix', default=False)
    depends_on('mpi')
        
    depends_on('netcdf-c@4.4.0 +shared +mpi')
    depends_on('hdf5@1.8.21 +fortran +shared +mpi')
    depends_on('pumi +int64 +zoltan +fortran', when='-simmetrix')
    depends_on('pumi +int64 simmodsuite=kernels +zoltan +fortran', when='+simmetrix')
    depends_on('zoltan@3.83 +parmetis+int64')

    def build_args(self, spec, prefix):                                                                               
        args=[]                                                                                                                                                                                                                         
        mpi_id = spec['mpi'].name + spec['mpi'].version.string                                                        
        args.append('mpiLib=' + mpi_id)                                                                               
        args.append('cc=mpicc')                                                                                       
        args.append('cxx=mpicxx')                                                                                                                                                                                                   
        if '+simModSuite' in spec:     
            args.append('simModSuite=yes')                                                                            
        return args                                                                                                   
                                                                                                                  
    def build(self, spec, prefix):                                                                                    
        args = self.build_args(spec, prefix)                                                                          
        scons(*args)                                                                                                  
    

    def install(self,spec,prefix):                                                                                  
        pass
        #make()
    