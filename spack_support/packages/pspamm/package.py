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
#     spack install pspamm
#
# You can edit this file again by typing:
#
#     spack edit pspamm
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *


class Pspamm(Package):
    homepage = "https://github.com/peterwauligmann/PSpaMM/blob/master/README.md"
    version('develop',
            git='https://github.com/peterwauligmann/PSpaMM',
            branch='master')

    variant('numpy', default=False, description="installs numpy")
    variant('scipy', default=False, description="installs scipy")

    depends_on('py-numpy', when='+numpy')
    depends_on('py-scipy', when='+scipy')

    def install(self, spec, prefix):
        install_tree('.', prefix)

    def setup_run_environment(self, env):
        env.prepend_path('PYTHONPATH', self.spec.prefix)
