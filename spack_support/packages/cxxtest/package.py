# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *

class Cxxtest(Package):
    """CxxTest is a unit testing framework for C++."""


    homepage = "https://http://cxxtest.com"
    version('develop',
            git='https://github.com/CxxTest/cxxtest',
            branch='master')

    def install(self, spec, prefix):
        install_tree('bin', prefix.bin)
        install_tree('cxxtest', prefix.include.cxxtest)
        install_tree('python', prefix.python)
        install_tree('doc', prefix.doc)
