# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""python -m seissol.trace MODEL.py [-o OUTPUT] [--model NAME] [--check]

Writes the core file of the model in MODEL.py (a function decorated with @sd.model), by default
MODEL.sderiv next to it, which SeisSol reads as python:MODEL.py. The model is imported and traced
in its own directory, so that a relative path in it is relative to it; the files it reads there
are recorded in the core file with their hashes, and SeisSol stops where one of them, or the
model, has changed since -- with the command that writes the file again.
"""

import argparse
import contextlib
import importlib.util
import os
import sys

from seissol.sderiv import Model, TracingError, recording


@contextlib.contextmanager
def _in_directory(directory):
    previous = os.getcwd()
    os.chdir(directory)
    sys.path.insert(0, directory)
    try:
        yield
    finally:
        sys.path.remove(directory)
        os.chdir(previous)


def _models(module):
    return {
        name: value for name, value in vars(module).items() if isinstance(value, Model)
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="python -m seissol.trace",
        description="Writes the core file of sderiv of a model in Python, for SeisSol.",
    )
    parser.add_argument("source", help="the file of the model, MODEL.py")
    parser.add_argument(
        "-o", "--output", help="the core file (by default MODEL.sderiv)"
    )
    parser.add_argument("--model", help="the model of the file, where it has several")
    parser.add_argument(
        "--check",
        action="store_true",
        help="load the core file with the module sderiv of a SeisSol build, where there is one",
    )
    args = parser.parse_args(argv)

    source = os.path.abspath(args.source)
    output = (
        os.path.abspath(args.output)
        if args.output
        else os.path.splitext(source)[0] + ".sderiv"
    )
    directory = os.path.dirname(output)
    command = (
        f"python -m seissol.trace {os.path.relpath(source, directory)} "
        f"-o {os.path.relpath(output, directory)}"
    )
    if args.model:
        command += f" --model {args.model}"

    try:
        with recording() as files, _in_directory(os.path.dirname(source)):
            name = "seissol_model_" + os.path.splitext(os.path.basename(source))[0]
            spec = importlib.util.spec_from_file_location(name, source)
            if spec is None:
                parser.error(f"{args.source} is no Python file")
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            models = _models(module)
            if args.model:
                if args.model not in models:
                    parser.error(
                        f"{args.source} has no model {args.model}; it has "
                        f"{', '.join(models) or 'none'}"
                    )
                chosen = models[args.model]
            elif len(models) == 1:
                chosen = next(iter(models.values()))
            else:
                parser.error(
                    f"{args.source} has {'no model' if not models else 'several models'}"
                    + (
                        f" ({', '.join(models)}): choose one with --model"
                        if models
                        else ""
                    )
                    + "; a model is a function decorated with @sd.model"
                )
            text = chosen.core(output, command=command, depends=files)
    except TracingError as error:
        print(f"seissol.trace: {error}", file=sys.stderr)
        return 1

    with open(output, "w", encoding="utf-8") as file:
        file.write(text)
    print(f"wrote {output}")

    if args.check:
        try:
            import sderiv
        except ImportError:
            print(
                "not checked: no module sderiv (a SeisSol build with -DSDERIV_PYTHON=ON)"
            )
            return 0
        program = sderiv.compile(text, file=output)
        print(
            f"program {program.fingerprint}: {', '.join(program.inputs)} -> "
            f"{', '.join(program.outputs)}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
