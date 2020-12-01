# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *
import os

class SeissolEnv(BundlePackage):
    """Seissol - A scientific software for the numerical simulation of seismic wave phenomena and earthquake dynamics.
    This package only provides all necessary libs for seissol installation.
    """

    homepage = "http://www.seissol.org"
    version('develop',
            git='https://github.com/SeisSol/SeisSol.git',
            branch='master')

    maintainers = ['ravil-mobile']
    
    variant('mpi', default=True, description="installs an MPI implementation")
    variant('asagi', default=True, description="installs asagi for material input")
    variant('extra_blas', default='none', description='installs an extra blas implementation',
            values=('mkl', 'openblas', 'blis', 'none'), 
            multi=True)
    variant('python', default=False, description="installs python, pip, numpy and scipy")
    variant('building_tools', default=False, description="installs scons and cmake")
    variant('x86', default=True, description="installs extra packages for x86 platform")
    

    depends_on('mpi', when="+mpi")
    depends_on('parmetis', when="+mpi")
    depends_on('metis +int64', when="+mpi")
    depends_on('libxsmm@1.15 +generator', when="+x86")

    depends_on('hdf5@1.8.21 +fortran +shared ~mpi', when="~mpi")
    depends_on('hdf5@1.8.21 +fortran +shared +mpi', when="+mpi")

    depends_on('netcdf-c@4.4.0 +shared ~mpi', when="~mpi")
    depends_on('netcdf-c@4.4.0 +shared +mpi', when="+mpi")

    depends_on('asagi ~mpi ~mpi3', when="+asagi ~mpi")
    depends_on('asagi +mpi +mpi3', when="+asagi +mpi")
    
    depends_on('intel-mkl threads=none', when="extra_blas=mkl")
    depends_on('openblas threads=none', when="extra_blas=openblas")
    depends_on('blis threads=none', when="extra_blas=blis")

    depends_on('memkind', when="+x86")
    depends_on('pspamm')
    depends_on('impalajit')
    depends_on('yaml-cpp@0.6.2')
    depends_on('cxxtest')
    

    
    depends_on('py-numpy', when='+python')
    depends_on('py-scipy', when='+python')
    depends_on('py-matplotlib', when='+python')
    depends_on('py-pip', when='+python')
    depends_on('py-pyopenssl', when='+python')
    depends_on('python@3.6.0', when='+python')
    

    depends_on('cmake@3.12.0:3.16.2', when='+building_tools')
    depends_on('scons@3.0.1:3.1.2', when='+building_tools')

    def setup_run_environment(self, env):
        
        roots = []; bins = []; libs = []; includes = []; pkgconfigs = []; pythonpath = []
        for child_spec in self.spec.dependencies():
            roots.append(child_spec.prefix if os.path.isdir(child_spec.prefix) else None)
            bins.append(child_spec.prefix.bin if os.path.isdir(child_spec.prefix.bin) else None)
            libs.append(child_spec.prefix.lib if os.path.isdir(child_spec.prefix.lib) else None)
            includes.append(child_spec.prefix.include if os.path.isdir(child_spec.prefix.include) else None)

            # one has to walk from the current root down in order to find pkgconfig folder
            # The reason is that some people include "pkgconfig" into "lib" but some put it into "share"
            # The second reason is to find all 'site-packages' and add them to PYTHONPATH
            for path, dirs, files in os.walk(child_spec.prefix):
                if "site-packages" in dirs:
                    pythonpath.append(os.path.join(path, "site-packages"))

                for file in files:
                    if file.endswith(".pc"):
                        pkgconfigs.append(path)
                        break

        env.prepend_path('CMAKE_PREFIX_PATH', ":".join(filter(None, roots)))
        env.prepend_path('PKG_CONFIG_PATH', ":".join(filter(None, pkgconfigs)))
        
        env.prepend_path('PATH', ":".join(filter(None, bins)))
        env.prepend_path('LD_LIBRARY_PATH', ":".join(filter(None, libs)))
        env.prepend_path('LIBRARY_PATH', ":".join(filter(None, libs)))
        
        env.prepend_path('CPATH', ":".join(filter(None, includes)))
        env.prepend_path('CPPPATH', ":".join(filter(None, includes)))
        env.prepend_path('C_INCLUDE_PATH', ":".join(filter(None, includes)))
        env.prepend_path('CPLUS_INCLUDE_PATH', ":".join(filter(None, includes)))

        env.prepend_path('PYTHONPATH', ":".join(filter(None, pythonpath)))
        
        # add pspamm while loading seissol-env
        env.prepend_path('PYTHONPATH', self.spec.dependencies_dict()['pspamm'].spec.prefix.pspamm)
        env.prepend_path('PATH', self.spec.dependencies_dict()['pspamm'].spec.prefix.pspamm)
