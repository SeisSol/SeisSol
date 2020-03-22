# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
import sys
import os
import glob

class SeissolUtils(Package):
    """Seissol - A scientific software for the numerical simulation of seismic wave phenomena 
    and earthquake dynamics. This package provides all necessary apps, libs, scripts 
    as well as simulation cases to work with SeisSol.
    """

    homepage = "http://www.seissol.org"
    version('develop',
            git='https://github.com/SeisSol/SeisSol.git',
            branch='ravil/spack', submodules=True)

    maintainers = ['ravil-mobile']

    variant('benchmarks', 
            default=False, 
            description="fetches benchmarks. Be sure that you "
                        "have access to the SeisSol LRZ-gitlab repo.")

    variant('paraview', default=False, description="installs Paraview for visualization")
    variant('building_tools', default=False, description="installs scons")

    resource(name='cookbook', 
             git='https://github.com/daisy20170101/SeisSol_Cookbook',
             placement='cookbook')

    resource(name='benchmarks', 
             git='https://gitlab.lrz.de/seissol.org/benchmarks',
             when='+benchmarks', 
             placement='benchmarks')

    depends_on("hdf5@1.8.21 +fortran +shared +mpi")
    depends_on('netcdf-c@4.4.0 +shared +mpi')
    depends_on('pumgen')
    depends_on('glm@0.9.7.1')
    depends_on('proj@4.9.2')

    depends_on("paraview", when="+paraview")    
    depends_on('scons@3.0.1:3.1.2', when='+building_tools')
    
    def install(self, spec, prefix):
        install_tree(self.stage.source_path.cookbook, prefix.cookbook)

        source_directory = self.stage.source_path
        # gmsh2gambit
        build_gmsh2gambit_dir = join_path(source_directory, 'preprocessing/meshing/gmsh2gambit')
        with working_dir(build_gmsh2gambit_dir, create=False):
            scons()
            files = glob.glob(join_path(build_gmsh2gambit_dir, 'build/bin', '*'))
            for file in files:
                install(file, prefix.bin)

        """
        #TODO: find a way how to install gambit2netcdf
        build_gambit2netcdf_dir = join_path(source_directory, 'preprocessing/meshing/gambit2netcdf')
        with working_dir(build_gambit2netcdf_dir, create=False):
            scons()
            files = glob.glob(join_path(build_gambit2netcdf_dir, 'build/bin', '*'))
            for file in files:
                install(file, prefix)
        """
        
        # cubes
        build_cube_dir = join_path(source_directory, 'preprocessing/meshing/cube_c')
        with working_dir(build_cube_dir, create=False):
            scons()
            files = glob.glob(join_path(build_cube_dir, 'build/bin', '*'))
            for file in files:
                install(file, prefix.bin)


        # rconv for point sources
        build_rconv_dir = join_path(source_directory, 'preprocessing/science/rconv')
        with working_dir(build_rconv_dir, create=False):
            scons('compiler={}'.format(*self.compiler.cc_names))
            files = glob.glob(join_path(build_rconv_dir, 'build/bin', '*'))
            for file in files:
                install(file, prefix.bin)
