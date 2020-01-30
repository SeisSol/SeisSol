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

    # FIXME: Add a proper url for your package's homepage here.
    homepage = "http://www.seissol.org"
    version('develop',
            git='https://github.com/SeisSol/SeisSol.git',
            branch='master')

    #version('master',  branch='master')
    maintainers = ['ravil-mobile']
    
    variant('mpi', default=True, description="use inter-node computing")
    variant('openmp', default=True, description="use intra-node computing")
    variant('asagi', default=True, description="use asagi for material input")


    depends_on('mpi', when="+mpi")

    #depends_on('parmetis ^metis+int64', when="+mpi")
    #depends_on('metis +int64+shared', when='+mpi')
    #depends_on('parmetis +shared', when='+mpi')

    depends_on('parmetis', when="+mpi")
    depends_on('metis+int64', when="~mpi")

    depends_on('libxsmm +generator')
    depends_on('memkind')

    depends_on('hdf5@1.8.21 +fortran +shared ~mpi', when="~mpi")
    depends_on('hdf5@1.8.21 +fortran +shared +mpi', when="+mpi")

    depends_on('netcdf-c@4.6.1 +shared ~mpi', when="~mpi")
    depends_on('netcdf-c@4.6.1 +shared +mpi', when="+mpi")

    depends_on('asagi ~mpi ~mpi3', when="+asagi ~mpi")
    depends_on('asagi +mpi +mpi3', when="+asagi +mpi")


    # instsall cxx_test manually
    # spack install seissol-env +mpi %gcc@8.3.0 ^openmpi@3.1.5 ^metis+int64
    # spack install seissol-env +mpi %intel@17.0.2 ^intel-mpi@2018.0.128
    # spack install seissol-env +mpi %intel@17.0.2 ^intel-mpi@2018.0.128 ^metis+int64
    #with working_dir("build", create=True):
    #   cmake("..", *std_cmake_args)
    #   make()