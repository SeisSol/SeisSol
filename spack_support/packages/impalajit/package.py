# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *

class Impalajit(CMakePackage):
    """A lightweight JIT compiler for flexible data access in simulation applications."""

    homepage = "https://github.com/uphoffc/ImpalaJIT/blob/master/README.md"
    version('develop',
            git='https://github.com/uphoffc/ImpalaJIT.git',
            branch='master',
            commit='0b2a2f503ab16')

    variant('static', default=True, description="compile as a static lib")
    
    # Tests are not working. They require the library being preinstalled
    # The option has been postponed until the dev. fixes the bug
    #variant('with_tests', default=True, description="compile and run tests")
    

    depends_on('cmake', type='build')
    depends_on('pkg-config', type='build')

    def cmake_args(self):
        args = []


        if '+static' in self.spec:
            args.append('-DSTATIC_LIB=ON')
            args.append('-DSHARED_LIB=OFF')
        else:
            args.append('-DSTATIC_LIB=OFF')
            args.append('-DSHARED_LIB=ON')

        # The option is not suported. See above
        #if '+with_tests' in self.spec:
        #    args.append('-DTESTS=ON')

        if self.compiler != "intel":
            args.append('-DINTEL_COMPILER=OFF')

        return args
