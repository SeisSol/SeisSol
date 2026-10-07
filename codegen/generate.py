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
import json
import os
import re
import sys

import kernels.arch
import kernels.configboundary
import kernels.dynamic_rupture
import kernels.general
import kernels.memlayout
import kernels.nodalbc
import kernels.output
import kernels.plasticity
import kernels.point
import kernels.quantities
import kernels.surface_displacement
import kernels.vtkproject
import yateto
from yateto import (
    DeviceArchDefinition,
    Generator,
    GlobalRoutineCache,
    HostArchDefinition,
    NamespacedGenerator,
    deriveArchitecture,
    fixArchitectureGlobal,
    gemm_configuration,
)
from yateto.ast.cost import BoundingBoxCostEstimator, FusedGemmsBoundingBoxCostEstimator
from yateto.gemm_configuration import GeneratorCollection
from yateto.metagen import MetaGenerator


def load_configs(cmdLineArgs):
    """The configurations to generate, in the order of their ids.

    Each one names the C++ type it is generated for (`key`) and the arguments it differs in from
    the command line. Without --configs, it is the one configuration that the command line
    describes, generated for seissol::Config0.
    """
    if cmdLineArgs.configs is not None:
        with open(cmdLineArgs.configs) as file:
            return json.load(file)
    return [
        {
            "key": "seissol::Config0",
            "equations": cmdLineArgs.equations,
            "order": cmdLineArgs.order,
            "precision": cmdLineArgs.precision,
            "numMechanisms": cmdLineArgs.numMechanisms,
            "multipleSimulations": cmdLineArgs.multipleSimulations,
            "drQuadRule": cmdLineArgs.drQuadRule,
            "solver": cmdLineArgs.solver,
        }
    ]


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
    cmdLineParser.add_argument(
        "--precision", type=str, choices=["s", "d", "f32", "f64"]
    )
    cmdLineParser.add_argument("--numMechanisms", type=int)
    cmdLineParser.add_argument("--vectorsize", default=0, type=int)
    cmdLineParser.add_argument("--alignment", default=0, type=int)
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
    cmdLineParser.add_argument(
        "--solver", type=str, choices=["linearck", "linearckanelastic", "stp"]
    )

    # "dry run" parameter for use directly in CMake (before building)
    cmdLineParser.add_argument(
        "--mode", type=str, choices=["collect", "codegen"], default="codegen"
    )
    cmdLineParser.add_argument("--codegen_target", type=str, default="__all__")
    # the configurations to generate, in the order of their ids; without it, the one that the
    # arguments above describe, under the key seissol::Config0
    cmdLineParser.add_argument("--configs", type=str, default=None)
    # the frameworks to write the unit tests of the kernels for, comma-separated, or none
    cmdLineParser.add_argument("--unit_tests", type=str, default="doctest")

    cmdLineParser.set_defaults(enable_premultiply_flux=False)
    cmdLineArgs = cmdLineParser.parse_args()

    # derive the compute platform
    gpu_platforms = ["cuda", "hip", "hipsycl", "acpp", "oneapi"]
    targets = ["gpu", "cpu"] if cmdLineArgs.device_backend in gpu_platforms else ["cpu"]

    if cmdLineArgs.vectorsize == 0:
        cmdLineArgs.vectorsize = None

    configs = load_configs(cmdLineArgs)
    # the arguments of each configuration: those of the command line, overridden by the fields of
    # the configuration
    configArgs = [
        argparse.Namespace(
            **{
                **vars(cmdLineArgs),
                **{key: value for key, value in config.items() if key != "key"},
            }
        )
        for config in configs
    ]

    def deriveWith(args, vectorsize):
        host = HostArchDefinition(args.host_arch, args.precision, vectorsize, None)
        device = None

        if args.device_backend != "none":
            device = DeviceArchDefinition(
                args.device_arch,
                args.device_vendor,
                args.device_backend,
                args.precision,
                vectorsize,
            )

        return deriveArchitecture(host, device), host, device

    def deriveArchitectureOf(args):
        arch, host_arch, device_arch = deriveWith(args, args.vectorsize)

        # The simulation index is the leading dimension of every fused tensor, and a
        # leading dimension is padded to the vector size. Padded simulation lanes
        # hold values nothing computes, and the hand-written parts of SeisSol index
        # the fused tensors with NumSimulations as the stride, so they would read
        # that padding as data. Narrow the vector size to the largest one the fused
        # simulations fill instead -- 32 B for eight single precision simulations on
        # a 64 B machine. The alignment a buffer starts on is a separate number and
        # keeps the architecture's value, which is why the two are derived apart.
        if args.multipleSimulations > 1:
            fusedBytes = args.multipleSimulations * arch.bytesPerReal
            vectorsize = arch.alignment
            while fusedBytes % vectorsize != 0:
                vectorsize //= 2
            if vectorsize != arch.alignment:
                print(
                    f"Reducing the vector size from {arch.alignment} B to "
                    f"{vectorsize} B, so that the {args.multipleSimulations} "
                    f"fused simulations are not padded.",
                    file=sys.stderr,
                )
                args.vectorsize = vectorsize
                arch, host_arch, device_arch = deriveWith(args, vectorsize)
        return arch, host_arch, device_arch

    archs = [deriveArchitectureOf(args) for args in configArgs]
    arch, host_arch, device_arch = archs[0]

    # One alignment and vector size hold for all configurations: the hand-written parts of SeisSol
    # align and pad their buffers by them.
    memoryCharacteristics = {
        (
            kernels.arch.cacheline(configArch),
            args.vectorsize or kernels.arch.vector_size(configArch),
        )
        for (configArch, _, _), args in zip(archs, configArgs)
    }
    if len(memoryCharacteristics) > 1:
        raise RuntimeError(
            "The configurations would need different alignments or vector sizes "
            f"(alignment, vector size in bytes: {sorted(memoryCharacteristics)}). "
            "Build them into executables of their own."
        )

    fixArchitectureGlobal(arch)

    os.makedirs(cmdLineArgs.outputDir, exist_ok=True)
    kernels.arch.emit_header(
        arch,
        cmdLineArgs.outputDir,
        override_alignment=cmdLineArgs.alignment,
        override_vectorsize=configArgs[0].vectorsize or 0,
    )

    # pick up the gemm tools defined by the user
    gemm_tool_list = re.split(r"[,;]", cmdLineArgs.gemm_tools.replace(" ", ""))

    def gemmToolsFor(arch):
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
        return GeneratorCollection(gemm_generators)

    cost_estimators = BoundingBoxCostEstimator
    custom_routine_generators = {}

    isOldGpuInterface = True

    if "gpu" in targets:
        device_codegen = re.split(r"[,;]", cmdLineArgs.device_codegen.replace(" ", ""))

        if "gemmforge-chainforge" in device_codegen and cmdLineArgs.device_backend in [
            "cuda",
            "hip",
        ]:
            chainforge_spec = importlib.util.find_spec("chainforge")
            if chainforge_spec is not None:
                cost_estimators = FusedGemmsBoundingBoxCostEstimator
            else:
                raise ModuleNotFoundError(
                    "Could not find chainforge. You can install it from github.com/seissol/chainforge ."
                )

        if "tensorforge" in device_codegen:
            import tensorforge

            isOldGpuInterface = False

            if tensorforge.use_fusedgemm_cost():
                cost_estimators = FusedGemmsBoundingBoxCostEstimator

            custom_routine_generators["gpu"] = tensorforge.get_routine_generator(yateto)

    subfolders = []

    routine_cache = GlobalRoutineCache()

    gemmTools = [gemmToolsFor(configArch) for configArch, _, _ in archs]

    # The code of the equation is named by the key of its configuration:
    # seissol::kernel::X<Cfg> is the kernel X of the configuration Cfg,
    # and seissol::Pool<Cfg> the pool it binds. runtime.h reaches the
    # kernels as well, by the variant of their configuration, with operands as
    # views where they are to take them so (see kernels.common.cold_kernel_attrs).
    metagen = MetaGenerator(["typename"])

    # Tensors that the code names whatever the configuration, but that only
    # some configurations have: in the others, the key names no tensor
    # (`void`), which kernels::size and kernels::familySize count as empty.
    optionalTensors = [
        "canonicalI",
        "coupledCanonicalI",
        "E",
        "ET",
        "Iane",
        "Qane",
        "Qext",
        "W",
        "Zinv",
        "dQane",
        "dQext",
        "normalStress",
        "spaceTimePredictor",
        "w",
    ]

    unitTests = (
        [] if cmdLineArgs.unit_tests == "none" else cmdLineArgs.unit_tests.split(",")
    )

    def check_run_codegen(name):
        return cmdLineArgs.mode == "codegen" and cmdLineArgs.codegen_target in (
            "__all__",
            name,
        )

    # the kernels of the faces between configurations, per configuration
    boundaryPlans = kernels.configboundary.plans(configArgs)
    # the conversion casts between precisions, which only TensorForge generates on GPUs
    boundaryTargets = [
        target for target in targets if target == "cpu" or not isOldGpuInterface
    ]

    def generate_equation(subfolders, args, arch, gemmTools, key, boundaryPlan):
        order = args.order
        # the tensors of the configuration are laid out for its architecture
        fixArchitectureGlobal(arch)
        precision = "double" if args.precision in ["d", "f64"] else "single"
        fusedSuffix = (
            "-f" + str(args.multipleSimulations) if args.multipleSimulations > 1 else ""
        )

        if args.memLayout == "auto":
            # TODO(Lukas) Don't hardcode this
            env = {
                "precision": args.precision,
                "equations": args.equations,
                "order": order,
                "arch": args.host_arch,
                "device_arch": args.device_arch,
                "multipleSimulations": args.multipleSimulations,
                "targets": targets,
                "gemmgen": gemm_tool_list,
            }
            mem_layout = kernels.memlayout.guessMemoryLayout(env)
        else:
            mem_layout = kernels.memlayout.resolveMemoryLayout(args.memLayout, targets)

        cmdArgsDict = dict(vars(args))
        cmdArgsDict["memLayout"] = mem_layout

        equationsModuleName = f"kernels.equations.{args.equations}"

        equationsSpec = importlib.util.find_spec(equationsModuleName)
        if equationsSpec is None:
            raise RuntimeError("Could not find kernels for " + args.equations)

        # actually load the module
        equations = importlib.import_module(equationsModuleName)

        equation_class = equations.kernel_class(**cmdArgsDict)

        adg = equation_class(**cmdArgsDict)

        include_tensors = set()
        generator = Generator(arch)

        # Equation-specific kernels
        adg.addInit(generator)
        adg.addLocal(generator, targets)
        adg.addNeighbor(generator, targets)
        adg.addTime(generator, targets)
        adg.add_include_tensors(include_tensors)

        kernels.vtkproject.addKernels(
            generator,
            adg,
            args.PlasticityMethod,
            args.matricesDir,
            targets,
        )
        kernels.vtkproject.includeTensors(args.matricesDir, include_tensors)

        # Common kernels
        include_tensors.update(
            kernels.dynamic_rupture.addKernels(
                NamespacedGenerator(generator, namespace="dynamicRupture"),
                adg,
                args.matricesDir,
                args.drQuadRule,
                targets,
                isOldGpuInterface,
            )
        )

        kernels.plasticity.addKernels(
            generator,
            adg,
            args.matricesDir,
            args.PlasticityMethod,
            targets,
        )
        kernels.plasticity.includeTensors(
            args.matricesDir, adg, args.PlasticityMethod, include_tensors
        )

        kernels.nodalbc.addKernels(
            generator,
            adg,
            include_tensors,
            args.matricesDir,
            args,
            targets,
        )
        kernels.surface_displacement.addKernels(
            generator, adg, include_tensors, targets
        )
        kernels.point.addKernels(generator, adg)

        riemannMaterial = kernels.configboundary.RIEMANN_MATERIAL[args.equations]
        kernels.configboundary.add_kernels(
            generator,
            adg,
            args.matricesDir,
            boundaryPlan,
            riemannMaterial,
            precision,
            boundaryTargets,
        )

        outputDirName = f"equation-{adg.name()}-{order}-{precision}{fusedSuffix}"
        # configurations that differ in other respects, e.g. the solver, get one each
        if outputDirName in subfolders:
            outputDirName += f"-{args.solver}-m{args.numMechanisms}-{args.drQuadRule}"
        if outputDirName in subfolders:
            raise RuntimeError(
                f"Two configurations would be generated into {outputDirName}."
            )
        trueOutputDir = os.path.join(args.outputDir, outputDirName)
        if not os.path.exists(trueOutputDir):
            os.mkdir(trueOutputDir)

        subfolders += [outputDirName]

        kernels.quantities.emit_header(adg, trueOutputDir, key)
        kernels.configboundary.emit_header(
            adg, trueOutputDir, key, boundaryPlan, riemannMaterial, boundaryTargets
        )

        metagen.add_generator(
            [key],
            generator,
            name=re.sub(r"\W", "_", outputDirName),
            directory=outputDirName,
            gemm_cfg=gemmTools,
            cost_estimator=cost_estimators,
            include_tensors=include_tensors,
            routine_exporters=custom_routine_generators,
            routine_cache=routine_cache,
            unit_tests=unitTests,
        )

        return outputDirName

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

        generator = Generator(arch)

        kernels.general.addStiffnessTensor(generator)
        kernels.dynamic_rupture.addKernelsGeneral(
            NamespacedGenerator(generator, namespace="dynamicRupture")
        )

        if check_run_codegen("general"):
            generator.generate(
                outputDir=outputDir,
                namespace="seissol::general",
                gemm_cfg=gemmTools[0],
                cost_estimator=cost_estimators,
                include_tensors=kernels.general.includeMatrices(
                    cmdLineArgs.matricesDir
                ),
                routine_exporters=custom_routine_generators,
                routine_cache=routine_cache,
                unit_tests=unitTests,
            )

    def forward_files(filename):
        # Not every subfolder emits every file: the quantity layout, for one,
        # only exists for the equation.
        present = [
            folder
            for folder in subfolders
            if os.path.exists(os.path.join(cmdLineArgs.outputDir, folder, filename))
        ]
        with open(os.path.join(cmdLineArgs.outputDir, filename), "w") as file:
            file.writelines(["// IWYU pragma: begin_exports\n"])
            file.writelines(
                [f'#include "{os.path.join(folder, filename)}"\n' for folder in present]
            )
            file.writelines(["// IWYU pragma: end_exports\n"])

    equationFolders = [
        generate_equation(
            subfolders, args, configArch, tools, config["key"], boundaryPlan
        )
        for args, (configArch, _, _), tools, config, boundaryPlan in zip(
            configArgs, archs, gemmTools, configs, boundaryPlans
        )
    ]

    # Generate code (if we need to): the metagen generates the code of all configurations at once
    if any(check_run_codegen(folder) for folder in equationFolders):
        metagen.generate(
            cmdLineArgs.outputDir,
            namespace="seissol",
            includes=["Config.h"],
            declarationsTensors=optionalTensors,
            # only some configurations have them (see kernels.configboundary.plans)
            declarationsKernels=[
                "toCanonical",
                "fromCanonical",
                "fromCoupledCanonical",
                "gpu_toCanonical",
                "gpu_fromCanonical",
                "gpu_fromCoupledCanonical",
            ],
        )
    generate_general(subfolders)

    if cmdLineArgs.mode == "codegen":
        routine_cache.generate(cmdLineArgs.outputDir, "seissol")

        # init.h, kernel.h, pool.h and tensor.h are the metagen's, which name
        # the code of the equation by key; the code of general/, which belongs
        # to no configuration, is included from there.
        forward_files("quantities.h")
        forward_files("configboundary.h")

    if cmdLineArgs.mode == "collect":
        doctests = (
            [f"{Generator.DOCTEST_FILE_NAME}.cpp"] if "doctest" in unitTests else []
        )
        targets = {
            folder: {
                "kernels": [
                    os.path.join(folder, "init.cpp"),
                    os.path.join(folder, "kernel.cpp"),
                    os.path.join(folder, "pool.cpp"),
                    os.path.join(folder, "tensor.cpp"),
                ],
                "tests": [os.path.join(folder, name) for name in doctests],
                "headers": [
                    os.path.join(folder, "init.h"),
                    os.path.join(folder, "kernel.h"),
                    os.path.join(folder, "pool.h"),
                    os.path.join(folder, "tensor.h"),
                ],
            }
            for folder in subfolders
        }

        # The metagen knows what it adds: a translation unit per generator that
        # binds views to its kernels, and those of runtime.h next to them.
        for sources in metagen.sources().values():
            targets[os.path.dirname(sources[0])]["kernels"] = sources
        targets["runtime"] = {
            "kernels": metagen.shared_sources(),
            "tests": [],
            "headers": metagen.shared_headers(),
        }

        kernels.output.write_if_changed(
            os.path.join(cmdLineArgs.outputDir, "targets.json"), json.dumps(targets)
        )


if __name__ == "__main__":
    main()
