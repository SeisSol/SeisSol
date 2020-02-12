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
#     spack install seissol-utils
#
# You can edit this file again by typing:
#
#     spack edit seissol-utils
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *
import sys
import os
import glob

class SeissolUtils(Package):
    """Seissol - A scientific software for the numerical simulation of seismic wave phenomena and earthquake dynamics.
    This package provides all necessary apps, libs, scripts, simulation cases to work with SeisSol.
    """

    homepage = "http://www.seissol.org"
    version('develop',
            git='https://github.com/SeisSol/SeisSol.git',
            branch='ravil/spack', submodules=True)

    maintainers = ['ravil-mobile']

    variant('mpi', default=True, description="use inter-node computing")
    variant('paraview', default=False, description="install Paraview for visualization")
    variant("with_simmetrix", default=False, description="enable SIMMETRIX support for PUMI")
    variant('with_benchmarks', default=False, 
            description="fetches benchmarks. Be sure that you have access to SeisSol LRZ-gitlab repo.")

    resource(name='cookbook', git='https://github.com/daisy20170101/SeisSol_Cookbook',
             placement='cookbook')

    resource(name='benchmarks', git='https://gitlab.lrz.de/seissol.org/benchmarks',
             when='+with_benchmarks', placement='benchmarks')

    #depends_on("hdf5@1.8.21 +fortran +shared +mpi")
    #depends_on("paraview", when="+paraview")    
    #depends_on("pumi +mpi simmodsuite=none", when="~with_simmetrix") # add -O2 -g -Wall
    #depends_on("pumi +mpi simmodsuite=base", when="+with_simmetrix") # add -O2 -g -Wall
    #depends_on('scons@3.0.1:3.1.2', when='+building_tools')
    depends_on('netcdf-c@4.6.1 +shared ~mpi', when="~mpi")
    depends_on('netcdf-c@4.6.1 +shared +mpi', when="+mpi")
    
    depends_on('glm@0.9.7.1')
    depends_on('proj@4.9.2')


    def install(self, spec, prefix):
        source_directory = self.stage.source_path

        build_gmsh2gambit_dir = join_path(source_directory, 'preprocessing/meshing/gmsh2gambit')
        with working_dir(build_gmsh2gambit_dir, create=False):
            scons()
            files = glob.glob(join_path(build_gmsh2gambit_dir, 'build/bin', '*'))
            for file in files:
                install(file, prefix)

        """
        build_gambit2netcdf_dir = join_path(source_directory, 'preprocessing/meshing/gambit2netcdf')
        with working_dir(build_gambit2netcdf_dir, create=False):
            scons()
            files = glob.glob(join_path(build_gambit2netcdf_dir, 'build/bin', '*'))
            for file in files:
                install(file, prefix)
        """

        build_cube_dir = join_path(source_directory, 'preprocessing/meshing/cube_c')
        with working_dir(build_cube_dir, create=False):
            scons()
            files = glob.glob(join_path(build_cube_dir, 'build/bin', '*'))
            for file in files:
                install(file, prefix)


        build_rconv_dir = join_path(source_directory, 'preprocessing/science/rconv')
        with working_dir(build_rconv_dir, create=False):
            scons('compiler={}'.format(*self.compiler.cc_names))
            files = glob.glob(join_path(build_rconv_dir, 'build/bin', '*'))
            for file in files:
                install(file, prefix)


        # install PUMGen from https://github.com/SeisSol/PUMGen.git + Zoltan
        #make()
        #make('install')
