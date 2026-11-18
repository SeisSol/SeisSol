# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The core text of a trace (docs/sderiv-spec.rst): the inputs, the grids, the states, a
definition per operation, the next values of the states, the outputs and the checks -- and the
provenance in front, which SeisSol checks when it loads the file: the source and the files read
while tracing with their hashes, the command that writes the file again, the hash of its body,
and where every definition and output comes from in the source."""

from __future__ import annotations

import hashlib
import os
import re

RESERVED = {
    "meta",
    "in",
    "grid",
    "def",
    "state",
    "next",
    "out",
    "check",
    "inf",
    "nan",
    "let",
    "type",
    "for",
    "import",
    "if",
    "then",
    "else",
    "and",
    "or",
    "not",
}

_PLAIN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")

_OPERATORS = {"add": "+", "sub": "-", "mul": "*", "div": "/"}


def spell_name(name):
    """A name as the core format writes it: plain, or in backquotes."""
    if _PLAIN.match(name) and name not in RESERVED:
        return name
    if not name or "`" in name or "\n" in name or "\r" in name:
        raise ValueError(f"the core format cannot write the name {name!r}")
    return f"`{name}`"


def spell_number(value):
    """A number, as Python writes it shortest: also inf, nan and -0.0."""
    return repr(float(value))


def spell_string(text):
    if "\n" in text or "\r" in text:
        raise ValueError(f"a string of the core format is on one line: {text!r}")
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _spell_value(value):
    if isinstance(value, str):
        return spell_string(value)
    if isinstance(value, (list, tuple)):
        return "(" + ", ".join(_spell_value(v) for v in value) + ")"
    return spell_number(value)


def _relative(path, directory):
    """`path` relative to `directory`, where both are; a name such as "<stdin>" as it is."""
    if directory is None or not os.path.isabs(path):
        return path
    return os.path.relpath(path, directory)


def spell_grid(name, grid, directory):
    """The declaration of `grid` as `name`; its path relative to `directory` (absolute where
    that is None)."""
    attributes = [("kind", grid.kind)]
    if grid.path is not None:
        attributes.append(("path", _relative(grid.path, directory)))
    if grid.variable is not None:
        attributes.append(("variable", grid.variable))
    attributes.append(("components", tuple(grid.components)))
    if grid.coordinates is not None:
        attributes.append(("coordinates", tuple(grid.coordinates)))
    if grid.axes is not None:
        attributes.append(("axes", tuple(grid.axes)))
    if grid.knots is not None:
        attributes.append(("knots", tuple(tuple(axis) for axis in grid.knots)))
        attributes.append(("values", _nested(grid.values)))
    for key, value in (
        ("interpolation", grid.interpolation),
        ("boundary", grid.boundary),
    ):
        attributes.append((key, value if isinstance(value, str) else tuple(value)))
    if grid.time is not None:
        attributes.append(("time", grid.time))
    spelled = ", ".join(f"{key} = {_spell_value(value)}" for key, value in attributes)
    return f"grid {spell_name(name)} ({spelled})"


def _nested(values):
    """The values of an inline grid, [axis 0]...[axis d-1][component], nested as the core format
    nests them: the last axis outermost, down to the components of a node."""
    dimensions = values.ndim - 1
    index = [0] * dimensions

    def nest(axis):
        if axis < 0:
            return tuple(float(v) for v in values[tuple(index)])
        nodes = []
        for i in range(values.shape[axis]):
            index[axis] = i
            nodes.append(nest(axis - 1))
        return tuple(nodes)

    return nest(dimensions - 1)


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as file:
        for block in iter(lambda: file.read(1 << 20), b""):
            digest.update(block)
    return "sha256:" + digest.hexdigest()


def core_text(trace, directory=None, generator=None, command=None, depends=()):
    """The core text of `trace`, with its provenance where `directory` (that of the core file)
    is given: the source of the model, the files `depends` (absolute paths), `command`, the hash
    of the body, and the positions in the source."""
    graph = trace.graph
    roots = [node for _, node, _ in trace.outputs]
    roots += [node for _, _, node, _ in trace.states]
    roots += [node for node, _, _ in graph.checks]

    # the operations the roots read, each after those it reads
    order = []
    seen = set()
    for root in roots:
        stack = [(root, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                order.append(node)
                continue
            if node in seen or graph.nodes[node].op in ("const", "input", "state"):
                continue
            seen.add(node)
            stack.append((node, True))
            for arg in reversed(graph.nodes[node].args):
                stack.append((arg, False))

    # names: the inputs and the states by theirs, the grids and the definitions by fresh ones
    taken = set(trace.inputs) | {name for name, _, _, _ in trace.states}
    counter = iter(range(1 << 62))

    def fresh(prefix):
        while True:
            name = f"{prefix}{next(counter)}"
            if name not in taken:
                taken.add(name)
                return name

    grid_names = [fresh("g") for _ in graph.grids]
    counter = iter(range(1 << 62))
    names = {}
    for node in order:
        names[node] = fresh("t")

    def operand(node):
        n = graph.nodes[node]
        if n.op == "const":
            return spell_number(n.value)
        if n.op in ("input", "state"):
            return spell_name(n.name)
        return names[node]

    def expression(node):
        n = graph.nodes[node]
        args = [operand(a) for a in n.args]
        if n.op in _OPERATORS:
            return f"{args[0]} {_OPERATORS[n.op]} {args[1]}"
        if n.op == "neg":
            return f"-{args[0]}"
        if n.op == "sample":
            return f"{grid_names[n.grid]}.{spell_name(n.name)}({', '.join(args)})"
        return f"{n.op}({', '.join(args)})"

    body = []
    if trace.inputs:
        body.append("in " + ", ".join(spell_name(name) for name in trace.inputs))
    for name, declared in zip(grid_names, graph.grids):
        body.append(spell_grid(name, declared, directory))
    for name, initial, _, _ in trace.states:
        body.append(f"state {spell_name(name)} = {spell_number(initial)}")
    for node in order:
        body.append(f"def {names[node]} = {expression(node)}")
    for name, _, node, _ in trace.states:
        body.append(f"next {spell_name(name)} = {operand(node)}")
    for name, node, _ in trace.outputs:
        body.append(f"out {spell_name(name)} = {operand(node)}")
    for node, message, _ in graph.checks:
        body.append(f"check {operand(node)}, {spell_string(message)}")

    meta = []
    if generator is not None:
        meta.append(("generator", [generator]))
    if directory is not None:
        if trace.source is not None:
            meta.append(
                (
                    "source",
                    [_relative(trace.source, directory), sha256_of(trace.source)],
                )
            )
        for path in sorted(set(depends), key=lambda p: _relative(p, directory)):
            meta.append(("depends", [_relative(path, directory), sha256_of(path)]))
        if command is not None:
            meta.append(("command", [command]))
        text = "".join(line + "\n" for line in body)
        meta.append(
            ("body", ["sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()])
        )
    # where the definitions and the outputs come from
    for node in order:
        at = graph.positions.get(node)
        if at is not None:
            meta.append(
                ("line", [names[node], _relative(at[0], directory), at[1], at[2]])
            )
    for name, _, at in trace.outputs:
        if at is not None:
            meta.append(("line", [name, _relative(at[0], directory), at[1], at[2]]))

    lines = ["#!sderiv 1"]
    for key, values in meta:
        spelled = ", ".join(
            spell_string(v) if isinstance(v, str) else str(int(v)) for v in values
        )
        lines.append(f"meta {key} = {spelled}")
    return "\n".join(lines + body) + "\n"
