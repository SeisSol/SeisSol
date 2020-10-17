# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


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
        install_tree('.', prefix.pspamm)

    def setup_run_environment(self, env):
        env.prepend_path('PYTHONPATH', self.spec.prefix.pspamm)
