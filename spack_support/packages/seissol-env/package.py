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
#     spack install example
#
# You can edit this file again by typing:
#
#     spack edit example
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *


class SeissolEnv(BundlePackage):
    """A scientific software for the numerical simulation of seismic wave phenomena and earthquake dynamics.
    This package only provides all necessary libs for seissol installation.
    """

    homepage = "http://www.seissol.org"
    version('develop',
            git='https://github.com/SeisSol/SeisSol.git',
            branch='master')

    maintainers = ['ravil-mobile']
    
    variant('mpi', default=True, description="use inter-node computing")
    variant('openmp', default=True, description="use intra-node computing")
    variant('asagi', default=True, description="use asagi for material input")

    variant('extra_blas', default='openblas', description='add an extra blas implementation along with libxsmm',
            values=('mkl', 'openblas', 'blis'), 
            multi=False)


    depends_on('mpi', when="+mpi")

    depends_on('parmetis', when="+mpi")
    depends_on('metis +int64', when="+mpi")

    depends_on('libxsmm +generator')
    depends_on('memkind')

    depends_on('hdf5@1.8.21 +fortran +shared ~mpi', when="~mpi")
    depends_on('hdf5@1.8.21 +fortran +shared +mpi', when="+mpi")

    depends_on('netcdf-c@4.6.1 +shared ~mpi', when="~mpi")
    depends_on('netcdf-c@4.6.1 +shared +mpi', when="+mpi")

    depends_on('asagi ~mpi ~mpi3', when="+asagi ~mpi")
    depends_on('asagi +mpi +mpi3', when="+asagi +mpi")
    
    #depends_on('intel-mkl -threads', when='+extra_blas=mkl')
    #depends_on('openblas -threads', when='+extra_blas=openblas')
    #depends_on('blis -threads', when='+extra_blas=blis')

    depends_on('pspamm')
    depends_on('impalajit')
    depends_on('yaml-cpp@0.6.2')
    depends_on('cxxtest')
