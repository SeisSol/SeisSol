# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Evaluate the easi models of the verification cases without easi.

Only what the cases use is understood: the components !ConstantMap,
!AffineMap, !Any, !Switch, !GroupFilter and !AxisAlignedCuboidalDomainFilter,
in block style, with flow style for lists of numbers and for the parameter
lists of !Switch. Anything else is refused, rather than read as something it
is not.

A component maps the values it receives, the coordinates x, y and z at the
root, to its own outputs and hands them to the first of its components that
accepts them, or returns them if it has none; a filter accepts what lies in
its domain and passes it on unchanged. !Switch takes each of its parameter
lists from the component it assigns to it. The inputs of !AffineMap are taken
in lexicographic order of their names, as in easi.
"""

import math
import re
from pathlib import Path

_COMPONENTS = (
    "ConstantMap",
    "AffineMap",
    "Any",
    "Switch",
    "GroupFilter",
    "AxisAlignedCuboidalDomainFilter",
)
_KEY = re.compile(r"^(\[[^\]]*\]|[A-Za-z_]\w*)\s*:(?:\s+(.*))?$")
_SPECIAL = {".inf": math.inf, "+.inf": math.inf, "-.inf": -math.inf}


class Unsupported(ValueError):
    """The model uses something this module does not understand."""


class Tagged:
    """A node of the model with the tag of its component."""

    def __init__(self, tag, value):
        if tag not in _COMPONENTS:
            raise Unsupported(f"!{tag} is not supported")
        self.tag = tag
        self.value = value


def _scalar(text):
    text = text.strip()
    if text.startswith("[") and text.endswith("]"):
        return [_scalar(item) for item in text[1:-1].split(",") if item.strip()]
    if text in _SPECIAL:
        return _SPECIAL[text]
    try:
        return float(text)
    except ValueError:
        return text


def _lines(path):
    lines = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        text = raw.split("#", 1)[0].rstrip()
        if text.strip():
            lines.append([len(text) - len(text.lstrip()), text.strip()])
    return lines


def _node(lines, index, indent):
    """The node whose first line is lines[index], at least `indent` deep."""
    tag = None
    if lines[index][1].startswith("!") and " " not in lines[index][1]:
        tag = lines[index][1][1:]
        index += 1
        if index == len(lines) or lines[index][0] < indent:
            return Tagged(tag, None), index
    depth, text = lines[index]
    if text.startswith("- ") or text == "-":
        value, index = _sequence(lines, index, depth)
    else:
        value, index = _mapping(lines, index, depth)
    return (Tagged(tag, value) if tag is not None else value), index


def _sequence(lines, index, depth):
    items = []
    while index < len(lines) and lines[index][0] == depth:
        text = lines[index][1]
        if not text.startswith("-"):
            break
        # the item continues as if its text started two columns further in
        lines[index] = [depth + 2, text[1:].strip()]
        item, index = _node(lines, index, depth + 2)
        items.append(item)
    return items, index


def _mapping(lines, index, depth):
    mapping = {}
    while index < len(lines) and lines[index][0] == depth:
        match = _KEY.match(lines[index][1])
        if match is None:
            raise Unsupported(f"expected name: value, got {lines[index][1]!r}")
        key, rest = match.group(1), (match.group(2) or "").strip()
        if key.startswith("["):
            key = tuple(str(item) for item in _scalar(key))
        index += 1
        deeper = index < len(lines) and lines[index][0] > depth
        if rest.startswith("!"):
            tag = rest[1:]
            if deeper:
                value, index = _node(lines, index, depth + 1)
                mapping[key] = Tagged(tag, value)
            else:
                mapping[key] = Tagged(tag, None)
        elif rest:
            mapping[key] = _scalar(rest)
        elif deeper:
            mapping[key], index = _node(lines, index, depth + 1)
        else:
            mapping[key] = None
    return mapping, index


def load(path):
    """The model of an easi file."""
    lines = _lines(path)
    if not lines:
        raise Unsupported(f"{path} is empty")
    node, index = _node(lines, 0, 0)
    if index != len(lines):
        raise Unsupported(f"{path}: cannot read {lines[index][1]!r}")
    if not isinstance(node, Tagged):
        raise Unsupported(f"{path} has no component at its root")
    return node


def evaluate(node, values, group=None):
    """The outputs of a component for the input `values` and the material
    group `group`; None where the input lies outside its domain."""
    body = node.value or {}
    if node.tag == "ConstantMap":
        outputs = {name: float(value) for name, value in body["map"].items()}
    elif node.tag == "AffineMap":
        names = sorted(values)
        outputs = {}
        for name, row in body["matrix"].items():
            if len(row) != len(names):
                raise Unsupported(f"!AffineMap row {name} does not match {names}")
            outputs[name] = sum(c * values[n] for c, n in zip(row, names))
            outputs[name] += float(body.get("translation", {}).get(name, 0.0))
    elif node.tag == "GroupFilter":
        groups = body["groups"]
        groups = groups if isinstance(groups, list) else [groups]
        if group not in [int(item) for item in groups]:
            return None
        outputs = values
    elif node.tag == "AxisAlignedCuboidalDomainFilter":
        for name, (lower, upper) in body["limits"].items():
            if not lower <= values[name] <= upper:
                return None
        outputs = values
    elif node.tag == "Switch":
        outputs = {}
        for names, component in body.items():
            result = evaluate(component, values, group)
            if result is None:
                return None
            outputs.update({name: result[name] for name in names})
        return outputs
    else:
        outputs = values
    components = body.get("components")
    if components is None:
        return outputs
    for component in components if isinstance(components, list) else [components]:
        result = evaluate(component, outputs, group)
        if result is not None:
            return result
    return None


def evaluate_at(node, point, group=None):
    """The parameters of the model at the point (x, y, z)."""
    return evaluate(node, dict(zip("xyz", (float(value) for value in point))), group)
