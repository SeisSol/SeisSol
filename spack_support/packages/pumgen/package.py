# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *

class Pumgen(CMakePackage):
    homepage = "https://github.com/SeisSol/PUMGen/wiki/How-to-compile-PUMGen"
    version('develop',
            git='https://github.com/SeisSol/PUMGen.git',
            branch='master',
            submodules=True)

    maintainers = ['ravil-mobile']
    variant('with_simmetrix', default=False)
    variant('with_netcdf', default=False)

    depends_on('mpi')
        
    depends_on('netcdf-c +shared +mpi', when='+with_netcdf') # NOTE: only tested with 4.4.0 version
    depends_on('hdf5 +fortran +shared +hl +mpi') # NOTE: only tested with 1.8.21 version
    depends_on('pumi +int64 +zoltan -fortran', when='~with_simmetrix')
    depends_on('simmetrix-simmodsuite', when='+with_simmetrix')
    depends_on('pumi +int64 simmodsuite=kernels +zoltan -fortran', when='+with_simmetrix')
    depends_on('zoltan@3.83 +parmetis+int64 -fortran')

    def cmake_args(self):
        args = [
            self.define_from_variant('SIMMETRIX', 'with_simmetrix'),
            self.define_from_variant('NETCDF', 'with_netcdf')
        ]
        return args                                                                                                 

    def install(self, spec, prefix):
        self.cmake(spec, prefix)
        self.build(spec, prefix)
        install_tree(self.build_directory, prefix.bin)
