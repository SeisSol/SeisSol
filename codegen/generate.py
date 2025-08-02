#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff
# SPDX-FileContributor: Sebastian Wolf

import argparse
import importlib.util
import os
import re
import sys

from typing import Union

import kernels.dynamic_rupture
import kernels.general
import kernels.memlayout
import kernels.nodalbc
import kernels.plasticity
import kernels.point
import kernels.surface_displacement
import kernels.vtkproject
import kernels.config
import yateto
from yateto import (
    Generator,
    GlobalRoutineCache,
    NamespacedGenerator,
    gemm_configuration,
    deriveArchitecture,
    HostArchDefinition,
    DeviceArchDefinition,
    fixArchitectureGlobal,
)
from yateto.ast.cost import BoundingBoxCostEstimator, FusedGemmsBoundingBoxCostEstimator
from yateto.gemm_configuration import GeneratorCollection

from yateto.metagen import MetaGenerator


def main():

    cmdLineParser = argparse.ArgumentParser()
    cmdLineParser.add_argument("--equations")
    cmdLineParser.add_argument("--matricesDir")
    cmdLineParser.add_argument("--outputDir")
    cmdLineParser.add_argument("--host_arch")
    cmdLineParser.add_argument("--device_backend", default=None)
    cmdLineParser.add_argument("--device_arch", default=None)
    cmdLineParser.add_argument("--device_vendor", default=None)
    cmdLineParser.add_argument("--order", type=int)
    cmdLineParser.add_argument("--precision", type=str, choices=["s", "d"])
    cmdLineParser.add_argument("--numberOfMechanisms", type=int)
    cmdLineParser.add_argument("--vectorsize", default=None, type=Union[None, int])
    cmdLineParser.add_argument("--memLayout")
    cmdLineParser.add_argument("--multipleSimulations", type=int)
    cmdLineParser.add_argument("--PlasticityMethod")
    cmdLineParser.add_argument("--gemm_tools")
    cmdLineParser.add_argument("--device_codegen")
    cmdLineParser.add_argument("--drQuadRule")
    cmdLineParser.add_argument("--enable_premultiply_flux", action="store_true")
    cmdLineParser.add_argument(
        "--disable_premultiply_flux",
        dest="enable_premultiply_flux",
        action="store_false",
    )
    cmdLineParser.add_argument("--executable_libxsmm", default="")
    cmdLineParser.add_argument("--executable_pspamm", default="")
    cmdLineParser.set_defaults(enable_premultiply_flux=False)
    cmdLineArgs = cmdLineParser.parse_args()

    # derive the compute platform
    gpu_platforms = ["cuda", "hip", "hipsycl", "acpp", "oneapi"]
    targets = ["gpu", "cpu"] if cmdLineArgs.device_backend in gpu_platforms else ["cpu"]

    host_arch = HostArchDefinition(
        cmdLineArgs.host_arch, cmdLineArgs.precision, cmdLineArgs.vectorsize, None
    )
    device_arch = None

    if cmdLineArgs.device_backend != "none":
        device_arch = DeviceArchDefinition(
            cmdLineArgs.device_arch,
            cmdLineArgs.device_vendor,
            cmdLineArgs.device_backend,
            cmdLineArgs.precision,
            cmdLineArgs.vectorsize,
        )

    arch = deriveArchitecture(host_arch, device_arch)
    fixArchitectureGlobal(arch)

    # pick up the gemm tools defined by the user
    gemm_tool_list = re.split(r"[,;]", cmdLineArgs.gemm_tools.replace(" ", ""))
    gemm_generators = []

    for tool in gemm_tool_list:
        if hasattr(gemm_configuration, tool):
            specific_gemm_class = getattr(gemm_configuration, tool)
            # take executable arguments, but only if they are not empty
            if (
                specific_gemm_class is gemm_configuration.LIBXSMM
                and cmdLineArgs.executable_libxsmm != ""
            ):
                gemm_generators.append(
                    specific_gemm_class(arch, cmdLineArgs.executable_libxsmm)
                )
            elif (
                specific_gemm_class is gemm_configuration.PSpaMM
                and cmdLineArgs.executable_pspamm != ""
            ):
                gemm_generators.append(
                    specific_gemm_class(arch, cmdLineArgs.executable_pspamm)
                )
            else:
                gemm_generators.append(specific_gemm_class(arch))
        elif tool.strip().lower() == "tensorforge":
            pass  # TODO: remove (hence differently placed than "none")
        elif tool.strip().lower() != "none":
            print(f'Unknown GEMM tool "{tool}". Please refer to the documentation.')
            sys.exit("failure")

    cost_estimators = BoundingBoxCostEstimator
    custom_routine_generators = {}
    if "gpu" in targets:
        device_codegen = re.split(r"[,;]", cmdLineArgs.device_codegen.replace(" ", ""))

        if "gemmforge-chainforge" in device_codegen and cmdLineArgs.device_backend in [
            "cuda",
            "hip",
        ]:
            chainforge_spec = importlib.util.find_spec("chainforge")
            if chainforge_spec is not None:
                chainforge_spec.loader.load_module()
                cost_estimators = FusedGemmsBoundingBoxCostEstimator
            else:
                raise ModuleNotFoundError(
                    "Could not find chainforge. You can install it from github.com/seissol/chainforge ."
                )
        if "tensorforge" in device_codegen:
            import tensorforge

            custom_routine_generators["gpu"] = tensorforge.get_routine_generator(yateto)

    subfolders = []

    routine_cache = GlobalRoutineCache()

    gemmTools = GeneratorCollection(gemm_generators)

    metagen = MetaGenerator(["typename"])

    viscomode = "None"
    if cmdLineArgs.equations == "viscoelastic":
        viscomode = "QuantityExtension"
        equations = "Viscoelastic"
    elif cmdLineArgs.equations == "viscoelastic2":
        viscomode = "AnelasticTensor"
        equations = "Viscoelastic"
    else:
        viscomode = "None"
        equations = cmdLineArgs.equations[0].upper() + cmdLineArgs.equations[1:]

    quadrule = cmdLineArgs.drQuadRule[0].upper() + cmdLineArgs.drQuadRule[1:]

    configs = [
        {
            "order": cmdLineArgs.order,
            "mechanisms": cmdLineArgs.numberOfMechanisms,
            "equation": equations,
            "precision": "F64" if cmdLineArgs.host_arch[:1] == "d" else "F32",
            "viscomode": viscomode,
            "drquadrule": quadrule,
            "numsims": cmdLineArgs.multipleSimulations,
        }
    ]
    configsTemp = [
        {
            "order": cmdLineArgs.order,
            "mechanisms": cmdLineArgs.numberOfMechanisms,
            "equation": cmdLineArgs.equations,
            "precision": "F64" if cmdLineArgs.host_arch[:1] == "d" else "F32",
            "viscomode": viscomode,
            "drquadrule": cmdLineArgs.drQuadRule,
            "numsims": cmdLineArgs.multipleSimulations,
        }
    ]

    def generate_equation(subfolders, config):
        equationsSpec = importlib.util.find_spec(
            f"kernels.equations.{config['equation']}"
        )
        if equationsSpec is None:
            raise RuntimeError("Could not find kernels for " + config["equation"])
        equations = equationsSpec.loader.load_module()

        equation_class = equations.EQUATION_CLASS

        if cmdLineArgs.memLayout == "auto":
            # TODO(Lukas) Don't hardcode this
            env = {
                "equations": config["equation"],
                "order": config["order"],
                "arch": cmdLineArgs.host_arch,
                "device_arch": cmdLineArgs.device_arch,
                "multipleSimulations": config["numsims"],
                "targets": targets,
                "gemmgen": gemm_tool_list,
            }
            mem_layout = kernels.memlayout.guessMemoryLayout(env)
        elif not os.path.isabs(cmdLineArgs.memLayout):
            print(
                f"Using the pre-defined memory layout config file {cmdLineArgs.memLayout}"
            )
            script_dir = os.path.dirname(os.path.abspath(__file__))
            mem_layout = os.path.join(script_dir, "config", cmdLineArgs.memLayout)
        else:
            print(f"Using the memory layout config file {cmdLineArgs.memLayout}")
            mem_layout = cmdLineArgs.memLayout

        cmdArgsDict = vars(cmdLineArgs)
        cmdArgsDict["memLayout"] = mem_layout

        cmdArgsDict["drQuadRule"] = config["drquadrule"]
        cmdArgsDict["multipleSimulations"] = config["numsims"]
        cmdArgsDict["mechanisms"] = config["mechanisms"]

        adg = equation_class(**cmdArgsDict)

        include_tensors = set()
        generator = Generator(arch)

        # Equation-specific kernels
        adg.addInit(generator)
        adg.addLocal(generator, targets)
        adg.addNeighbor(generator, targets)
        adg.addTime(generator, targets)
        adg.add_include_tensors(include_tensors)

        kernels.vtkproject.addKernels(generator, adg, cmdLineArgs.matricesDir, targets)
        kernels.vtkproject.includeTensors(cmdLineArgs.matricesDir, include_tensors)

        # Common kernels
        include_tensors.update(
            kernels.dynamic_rupture.addKernels(
                NamespacedGenerator(generator, namespace="dynamicRupture"),
                adg,
                cmdLineArgs.matricesDir,
                cmdLineArgs.drQuadRule,
                targets,
            )
        )

        kernels.plasticity.addKernels(
            generator,
            adg,
            cmdLineArgs.matricesDir,
            cmdLineArgs.PlasticityMethod,
            targets,
        )
        kernels.nodalbc.addKernels(
            generator,
            adg,
            include_tensors,
            cmdLineArgs.matricesDir,
            cmdLineArgs,
            targets,
        )
        kernels.surface_displacement.addKernels(
            generator, adg, include_tensors, targets
        )
        kernels.point.addKernels(generator, adg)

        metagen.add_generator(
            ["Config"],
            generator,
            gemm_cfg=gemmTools,
            cost_estimator=cost_estimators,
            include_tensors=include_tensors,
            routine_exporters=custom_routine_generators,
            routine_cache=routine_cache,
        )

    def generate_general(subfolders):
        # we use always use double here,
        # since these kernels are only used in the initialization
        new_host_arch = HostArchDefinition(
            host_arch.archname, "d", host_arch.alignment, host_arch.prefetch
        )
        arch = deriveArchitecture(new_host_arch, None)
        fixArchitectureGlobal(arch)

        outputDir = os.path.join(cmdLineArgs.outputDir, "general")
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)

        subfolders += ["general"]

        # for now, enforce Eigen as a code generator here...
        # ...until we have a shared subroutine cache
        generator = Generator(arch)
        kernels.general.addStiffnessTensor(generator)
        generator.generate(
            outputDir=outputDir,
            namespace="seissol_general",
            gemm_cfg=gemmTools,
            cost_estimator=cost_estimators,
            include_tensors=kernels.general.includeMatrices(cmdLineArgs.matricesDir),
            routine_exporters=custom_routine_generators,
            routine_cache=routine_cache,
        )

    def forward_files(filename):
        with open(os.path.join(cmdLineArgs.outputDir, filename), "w") as file:
            file.writelines(["// IWYU pragma: begin_exports\n"])
            file.writelines(
                [
                    f'#include "{os.path.join(folder, filename)}"\n'
                    for folder in subfolders
                ]
            )
            file.writelines(["// IWYU pragma: end_exports\n"])

    for config in configsTemp:
        generate_equation(subfolders, config)
    generate_general(subfolders)

    with open(os.path.join(cmdLineArgs.outputDir, "..", "Config.h"), "w") as file:
        file.write(kernels.config.make_configfile(configs))

    with open(
        os.path.join(cmdLineArgs.outputDir, "..", "ConfigInclude.h"), "w"
    ) as file:
        file.write(kernels.config.make_configincludefile(configs))

    metagen.generate(
        os.path.join(cmdLineArgs.outputDir, "metagen"), "seissol", ["Config.h"]
    )
    subfolders += ["metagen"]

    routine_cache.generate(cmdLineArgs.outputDir, "seissol")

    forward_files("init.h")
    forward_files("kernel.h")
    forward_files("tensor.h")
    forward_files("init.cpp")
    forward_files("kernel.cpp")
    forward_files("tensor.cpp")
    forward_files("test-kernel.cpp")


if __name__ == "__main__":
    main()
