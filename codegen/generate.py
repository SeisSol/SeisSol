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

# The code is generated in steps, which a build runs as commands of their own, and so at the same
# time where it can. Every configuration is a step: it writes its code into its directory, and
# what the shared step needs to know of it -- what the metagen reports, and its routines, written
# already -- into EXCHANGE_FILE there. The shared step writes what the configurations share: the
# headers of the metagen, general/, the routines of all of them, and the headers that include
# those of every configuration. --mode collect lists the files of every step in steps.json.
STEP_ALL = "__all__"
STEP_SHARED = "__shared__"
EXCHANGE_FILE = "codegen.json"


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


def config_folders(configArgs):
    """The directory of each configuration, below the output directory."""
    folders = []
    for args in configArgs:
        precision = "double" if args.precision in ["d", "f64"] else "single"
        fusedSuffix = (
            "-f" + str(args.multipleSimulations) if args.multipleSimulations > 1 else ""
        )
        folder = f"equation-{args.equations}-{args.order}-{precision}{fusedSuffix}"
        # configurations that differ in other respects, e.g. the solver, get one each
        if folder in folders:
            folder += f"-{args.solver}-m{args.numMechanisms}-{args.drQuadRule}"
        if folder in folders:
            raise RuntimeError(f"Two configurations would be generated into {folder}.")
        folders.append(folder)
    return folders


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
    # what to generate (see STEPS): a configuration, by its directory; what the configurations share
    # (__shared__); or all of it, one after the other (__all__)
    cmdLineParser.add_argument("--codegen_target", type=str, default=STEP_ALL)
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
    # not by the steps of the configurations, which run at the same time
    if cmdLineArgs.mode == "collect" or cmdLineArgs.codegen_target in (
        STEP_ALL,
        STEP_SHARED,
    ):
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

    routine_cache = GlobalRoutineCache()

    gemmTools = [gemmToolsFor(configArch) for configArch, _, _ in archs]

    # The code of the equation is named by the key of its configuration:
    # seissol::kernel::X<Cfg> is the kernel X of the configuration Cfg,
    # and seissol::Pool<Cfg> the pool it binds. runtime.h reaches the
    # kernels as well, by the variant of their configuration, with operands as
    # views where they are to take them so (see kernels.common.cold_kernel_attrs).
    # Each configuration is generated in a step of its own (see STEP_ALL), and
    # the shared step writes the headers of the metagen from what they report:
    # the metagen needs no generator for that.
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

    # the kernels of the faces between configurations, per configuration
    boundaryPlans = kernels.configboundary.plans(configArgs)
    # the conversion casts between precisions, which only TensorForge generates on GPUs
    boundaryTargets = [
        target for target in targets if target == "cpu" or not isOldGpuInterface
    ]

    folders = config_folders(configArgs)
    for folder, config in zip(folders, configs):
        metagen.add_generator(
            [config["key"]],
            None,
            name=re.sub(r"\W", "_", folder),
            directory=folder,
            routine_cache=routine_cache,
            unit_tests=unitTests,
        )

    def generate_equation(args, arch, gemmTools, key, boundaryPlan, folder):
        order = args.order
        # the tensors of the configuration are laid out for its architecture
        fixArchitectureGlobal(arch)
        precision = "double" if args.precision in ["d", "f64"] else "single"

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

        trueOutputDir = os.path.join(args.outputDir, folder)
        os.makedirs(trueOutputDir, exist_ok=True)

        kernels.quantities.emit_header(adg, trueOutputDir, key)
        kernels.configboundary.emit_header(
            adg, trueOutputDir, key, boundaryPlan, riemannMaterial, boundaryTargets
        )

        return generator, include_tensors

    def generate_configuration(index):
        """The step of a configuration: its code, and EXCHANGE_FILE."""
        args, (configArch, _, _), tools, config, boundaryPlan, folder = list(
            zip(configArgs, archs, gemmTools, configs, boundaryPlans, folders)
        )[index]
        generator, include_tensors = generate_equation(
            args, configArch, tools, config["key"], boundaryPlan, folder
        )

        # the routines of this configuration alone
        cache = GlobalRoutineCache()
        single = MetaGenerator(["typename"])
        single.add_generator(
            [config["key"]],
            generator,
            name=re.sub(r"\W", "_", folder),
            directory=folder,
            gemm_cfg=tools,
            cost_estimator=cost_estimators,
            include_tensors=include_tensors,
            routine_exporters=custom_routine_generators,
            routine_cache=cache,
            unit_tests=unitTests,
        )
        summary = single.generate_single(0, cmdLineArgs.outputDir, namespace="seissol")
        with open(
            os.path.join(cmdLineArgs.outputDir, folder, EXCHANGE_FILE), "w"
        ) as file:
            # the directories relative to the output directory, so that the file is
            # the same wherever that is
            json.dump(
                {
                    "summary": summary,
                    "routines": cache.export(root=cmdLineArgs.outputDir),
                },
                file,
            )

    def generate_general():
        # we use always use double here,
        # since these kernels are only used in the initialization
        new_host_arch = HostArchDefinition(
            host_arch.archname, "d", host_arch.alignment, host_arch.prefetch
        )
        arch = deriveArchitecture(new_host_arch, None)
        fixArchitectureGlobal(arch)

        outputDir = os.path.join(cmdLineArgs.outputDir, "general")
        os.makedirs(outputDir, exist_ok=True)

        generator = Generator(arch)

        kernels.general.addStiffnessTensor(generator)
        kernels.dynamic_rupture.addKernelsGeneral(
            NamespacedGenerator(generator, namespace="dynamicRupture")
        )

        # written already, as those of the configurations are
        cache = GlobalRoutineCache()
        generator.generate(
            outputDir=outputDir,
            namespace="seissol::general",
            gemm_cfg=gemmTools[0],
            cost_estimator=cost_estimators,
            include_tensors=kernels.general.includeMatrices(cmdLineArgs.matricesDir),
            routine_exporters=custom_routine_generators,
            routine_cache=cache,
            unit_tests=unitTests,
        )
        return cache.export()

    def forward_files(filename):
        # Not every subfolder emits every file: the quantity layout, for one,
        # only exists for the equation.
        present = [
            folder
            for folder in folders + ["general"]
            if os.path.exists(os.path.join(cmdLineArgs.outputDir, folder, filename))
        ]
        with open(os.path.join(cmdLineArgs.outputDir, filename), "w") as file:
            file.writelines(["// IWYU pragma: begin_exports\n"])
            file.writelines(
                [f'#include "{os.path.join(folder, filename)}"\n' for folder in present]
            )
            file.writelines(["// IWYU pragma: end_exports\n"])

    def generate_shared():
        """The shared step, from EXCHANGE_FILE of every configuration."""
        exchanges = []
        for folder in folders:
            with open(
                os.path.join(cmdLineArgs.outputDir, folder, EXCHANGE_FILE)
            ) as file:
                exchanges.append(json.load(file))

        # the metagen generates the headers of all configurations at once
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
            precompiled=[exchange["summary"] for exchange in exchanges],
        )

        for exchange in exchanges:
            routine_cache.merge(exchange["routines"], root=cmdLineArgs.outputDir)
        routine_cache.merge(generate_general())
        routine_cache.generate(cmdLineArgs.outputDir, "seissol")

        # init.h, kernel.h, pool.h and tensor.h are the metagen's, which name
        # the code of the equation by key; the code of general/, which belongs
        # to no configuration, is included from there.
        forward_files("quantities.h")
        forward_files("configboundary.h")

    if cmdLineArgs.mode == "collect":
        # what Generator.generate writes for a generator in its directory
        def generator_files(folder):
            files = [
                os.path.join(folder, f"{name}.{extension}")
                for name in (
                    Generator.TENSORS_FILE_NAME,
                    Generator.INIT_FILE_NAME,
                    Generator.KERNELS_FILE_NAME,
                    Generator.POOL_FILE_NAME,
                )
                for extension in ("h", "cpp")
            ]
            if "doctest" in unitTests:
                files += [os.path.join(folder, f"{Generator.DOCTEST_FILE_NAME}.cpp")]
            if "cxxtest" in unitTests:
                files += [os.path.join(folder, f"{Generator.CXXTEST_FILE_NAME}.h")]
            # the header of its routines, which includes the one of all routines
            files += [os.path.join(folder, f"{Generator.ROUTINES_FILE_NAME}.h")]
            return files

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
            for folder in folders + ["general"]
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

        # The files of every step, and the files of other steps it reads. The
        # shared step writes into the directories of the configurations as
        # well: the units that bind views to their kernels for runtime.h, and
        # the headers of their routines.
        steps = [
            {
                "target": folder,
                "outputs": [
                    file
                    for file in generator_files(folder)
                    if not file.endswith(f"{Generator.ROUTINES_FILE_NAME}.h")
                ]
                + [
                    os.path.join(folder, "quantities.h"),
                    os.path.join(folder, "configboundary.h"),
                    os.path.join(folder, EXCHANGE_FILE),
                ],
                "inputs": [],
            }
            for folder in folders
        ]
        steps += [
            {
                "target": STEP_SHARED,
                "outputs": metagen.shared_headers()
                + metagen.shared_sources()
                + [
                    os.path.join(folder, f"{MetaGenerator.RUNTIME_NAME}.cpp")
                    for folder in folders
                ]
                + [
                    os.path.join(folder, f"{Generator.ROUTINES_FILE_NAME}.h")
                    for folder in folders
                ]
                + generator_files("general")
                + [
                    f"{Generator.ROUTINES_FILE_NAME}.h",
                    f"{Generator.ROUTINES_FILE_NAME}.cpp",
                    f"{Generator.GPULIKE_ROUTINES_FILE_NAME}.cpp",
                    "quantities.h",
                    "configboundary.h",
                ],
                "inputs": [os.path.join(folder, EXCHANGE_FILE) for folder in folders],
            }
        ]

        kernels.output.write_if_changed(
            os.path.join(cmdLineArgs.outputDir, "targets.json"), json.dumps(targets)
        )
        kernels.output.write_if_changed(
            os.path.join(cmdLineArgs.outputDir, "steps.json"), json.dumps(steps)
        )
        return

    target = cmdLineArgs.codegen_target
    if target not in folders + [STEP_SHARED, STEP_ALL]:
        raise RuntimeError(
            f'Unknown code generation target "{target}"; there are {STEP_ALL}, '
            f"{STEP_SHARED} and the directories of the configurations, "
            f"{', '.join(folders)}."
        )
    for index, folder in enumerate(folders):
        if target in (STEP_ALL, folder):
            generate_configuration(index)
    if target in (STEP_ALL, STEP_SHARED):
        generate_shared()


if __name__ == "__main__":
    main()
