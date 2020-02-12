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
#     spack install cxxtest
#
# You can edit this file again by typing:
#
#     spack edit cxxtest
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

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
