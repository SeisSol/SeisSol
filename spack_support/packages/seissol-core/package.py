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
#     spack install seissol-core
#
# You can edit this file again by typing:
#
#     spack edit seissol-core
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *


class SeissolCore(CMakePackage):
    """A scientific software for the numerical simulation of seismic wave phenomena and earthquake dynamics."""

    homepage = "http://www.seissol.org"
    version('develop',
            git='https://github.com/SeisSol/SeisSol.git',
            branch='master',
            submodules=True)

    maintainers = ['ravil-mobile']
    
    variant('gemm_tools', 
            default='LIBXSMM,PSpaMM', 
            description='Gemm tools for compute kernels',
            values=('LIBXSMM', 'PSpaMM', 'MKL', 'OpenBLAS', 'ACL_DEVICE_BLAS', 'ACL_DEVICE_BLAS'), 
            multi=False
    )
    variant('mpi', default=False)

    depends_on('seissol-env')
    depends_on('cmake', type='build')
    depends_on('py-numpy', type='run')

    def install(self, spec, prefix):
        # FIXME: Unknown build system
        make()
        make('install')
