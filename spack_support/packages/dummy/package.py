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
#     spack install dummy
#
# You can edit this file again by typing:
#
#     spack edit dummy
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *
import sys
import os
import glob

class Dummy(Package):
    homepage = "http://www.seissol.org"
    version('develop',
            git='https://github.com/SeisSol/SeisSol.git',
            branch='ravil/spack', submodules=True)

    variant('numpy', default=False, description="installes numpy if your platform doen't have numpy")
    variant('scipy', default=False, description="installes scipy if your platform doen't have scipy")
    variant('with_benchmarks', default=False, 
            description="fetches benchmarks. Be sure that you have access to SeisSol LRZ-gitlab repo.")
    variant('simmetrix', default=False)

    depends_on('py-numpy', when='+numpy')
    depends_on('py-scipy', when='+scipy')

    depends_on('mpi')
    #depends_on('hdf5@1.8.21 +fortran +shared +mpi')
    #depends_on('netcdf-c@4.4.0 +shared +mpi')
    #depends_on('glm@0.9.7.1')
    #depends_on('proj@4.9.2')
    depends_on('pumi +int64 ^zoltan@3.8.3', when='-simmetrix')
    depends_on('pumi +int64 simmodsuite=kernels ^zoltan@3.8.3', when='+simmetrix')
    
    #depends_on('impalajit')
    #depends_on('pspamm')
    #depends_on('asagi -mpi -mpi3 -numa')

    #resource(name='cookbook', git='https://github.com/daisy20170101/SeisSol_Cookbook',
    #         placement='cookbook')

    #resource(name='benchmarks', git='https://gitlab.lrz.de/seissol.org/benchmarks',
    #         when='+with_benchmarks', placement='benchmarks')

    #variant('extra_blas', default='none', description='add an extra blas implementation along with libxsmm',
    #        values=('mkl', 'openblas', 'blis', 'none'), 
    #        multi=True)

    #depends_on('openblas threads=none', when="extra_blas=openblas")
    #depends_on('blis threads=none', when="extra_blas=blis")
    #depends_on('intel-mkl threads=none', when="extra_blas=mkl")
    
    def install(self, spec, prefix):
        source_directory = self.stage.source_path
        build_rconv_dir = join_path(source_directory, 'preprocessing/science/rconv')
        with working_dir(build_rconv_dir, create=False):
            scons('compiler={}'.format(*self.compiler.cc_names))
            files = glob.glob(join_path(build_rconv_dir, 'build/bin', '*'))
            for file in files:
                install(file, prefix)


    def setup_run_environment(self, env):
        
        """
        roots = []; bins = []; libs = []; includes = []; pkgconfigs = []
        for child_spec in self.spec.dependencies():
            roots.append(child_spec.prefix if os.path.isdir(child_spec.prefix) else None)
            bins.append(child_spec.prefix.bin if os.path.isdir(child_spec.prefix.bin) else None)
            libs.append(child_spec.prefix.lib if os.path.isdir(child_spec.prefix.lib) else None)
            includes.append(child_spec.prefix.include if os.path.isdir(child_spec.prefix.include) else None)

            # one has to walk from the current root down in order to find pkgconfig folder
            # The reason is that some people include "pkgconfig" into "lib" but some put it into "share"
            for path, dirs, files in os.walk(child_spec.prefix):
                for file in files:
                    if file.endswith(".pc"):
                        pkgconfigs.append(path)
                        break


     
        env.prepend_path('MY_CMAKE_PREFIX_PATH', ":".join(filter(None, roots)) )
        env.prepend_path('MY_PATH', ":".join(filter(None, bins)) )
        env.prepend_path('MY_LIBS', ":".join(filter(None, libs)) )

        env.prepend_path('MY_CPATH', ":".join(filter(None, includes)) )
        #env.prepend_path('MY_CPPPATH', ":".join(filter(None, includes)) )
        #env.prepend_path('MY_C_INCLUDE_PATH', ":".join(filter(None, includes)) )
        #env.prepend_path('MY_CPLUS_INCLUDE_PATH', ":".join(filter(None, includes)) )
        env.prepend_path('MY_PKG_CONFIG_PATH', ":".join(filter(None, pkgconfigs)) )
        """
        pass


    #def install(self, spec, prefix):
    #    install_tree('.', prefix.bin)
    #    #raise	 