# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The conversion of a model: its branches rewritten so that a condition on a value of the point
runs both of them and merges what they leave behind (_branches.py), as AutoGraph converts a
function for TensorFlow.

An if becomes a call with its branches as functions of their own, which bind the names of the
model (nonlocal); a return, a break or a continue in a branch sets a flag, and the rest of the
block runs where it is not set; the operands of a conditional expression, and, or and a chained
comparison become thunks; not, min(), max(), any(), all(), assert and raise calls of their own.
The function is parsed from its source and compiled again at its lines, so that a trace names the
lines of the model, with its globals, its closure and its defaults; the functions defined in it
are converted with it. A construct the conversion cannot follow leaves the function as it is
written, with a warning.
"""

import __future__

import ast
import inspect
import sys
import types
import warnings

from . import _branches

_RT = "_sd__rt"
_RETVAL = _branches.RETURN
_RETURNED = "_sd__returned"
_PREFIX = "_sd__"

_DEFINITIONS = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
_LOOPS = (ast.For, ast.AsyncFor, ast.While)
_TRY = (ast.Try, ast.TryStar) if hasattr(ast, "TryStar") else (ast.Try,)

# the methods that change a container in place, which a branch saves it for
_MUTATORS = frozenset(
    {
        "append",
        "extend",
        "insert",
        "remove",
        "pop",
        "popitem",
        "clear",
        "update",
        "setdefault",
        "sort",
        "reverse",
        "add",
        "discard",
        "fill",
        "put",
        "itemset",
        "resize",
        "__setitem__",
        "__delitem__",
    }
)

_EXTREMA = ("min", "max", "any", "all")


class Unconvertible(Exception):
    """A construct the conversion cannot follow, and its node."""

    def __init__(self, message, node):
        super().__init__(message)
        self.node = node


# ---------------------------------------------------------------- what the source holds


def _blocks(stmt):
    """The blocks of statements of `stmt`."""
    if isinstance(stmt, (ast.If,) + _LOOPS):
        return [stmt.body, stmt.orelse]
    if isinstance(stmt, (ast.With, ast.AsyncWith)):
        return [stmt.body]
    if isinstance(stmt, _TRY):
        return [
            stmt.body,
            *(h.body for h in stmt.handlers),
            stmt.orelse,
            stmt.finalbody,
        ]
    if isinstance(stmt, ast.Match):
        return [case.body for case in stmt.cases]
    return []


def _statements(stmts, loops=True):
    """The statements of `stmts` and of their blocks in the same function; those in the body of a
    loop only where `loops` (the else of a loop is outside of it)."""
    for stmt in stmts:
        yield stmt
        if isinstance(stmt, _DEFINITIONS):
            continue
        for index, block in enumerate(_blocks(stmt)):
            if not loops and isinstance(stmt, _LOOPS) and index == 0:
                continue
            yield from _statements(block, loops)


def _has(kind, stmts, loops=True):
    return any(isinstance(stmt, kind) for stmt in _statements(stmts, loops))


def _leaves(stmts):
    """Whether `stmts` never go on after them: every way through them ends in a return, a raise,
    a break or a continue."""
    for stmt in stmts:
        if isinstance(stmt, (ast.Return, ast.Raise, ast.Break, ast.Continue)):
            return True
        if isinstance(stmt, ast.If) and _leaves(stmt.body) and _leaves(stmt.orelse):
            return True
        if isinstance(stmt, (ast.With, ast.AsyncWith)) and _leaves(stmt.body):
            return True
        if isinstance(stmt, _TRY):
            if _leaves(stmt.finalbody):
                return True
            if (_leaves(stmt.body) or _leaves(stmt.orelse)) and all(
                _leaves(handler.body) for handler in stmt.handlers
            ):
                return True
    return False


def _nodes(nodes, comprehensions=True):
    """The nodes of `nodes` and under them in the same scope: not in a function, a class or a
    lambda, nor in the targets of a comprehension, whose names are its own -- nor in a
    comprehension at all, but where `comprehensions`."""
    stack = list(nodes)
    while stack:
        node = stack.pop()
        yield node
        if isinstance(node, _DEFINITIONS + (ast.Lambda,)):
            continue
        if isinstance(
            node, (ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp)
        ):
            if not comprehensions:
                continue
        if isinstance(node, ast.comprehension):
            stack.append(node.iter)
            stack.extend(node.ifs)
            continue
        stack.extend(ast.iter_child_nodes(node))


def _bound(stmts):
    """The names `stmts` bind in their function (but those of an except, which it unbinds)."""
    names = set()
    for node in _nodes(stmts):
        if isinstance(node, ast.Name) and isinstance(node.ctx, (ast.Store, ast.Del)):
            names.add(node.id)
        elif isinstance(node, _DEFINITIONS):
            names.add(node.name)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            names.update((a.asname or a.name).split(".")[0] for a in node.names)
        elif isinstance(node, (ast.MatchAs, ast.MatchStar)) and node.name:
            names.add(node.name)
        elif isinstance(node, ast.MatchMapping) and node.rest:
            names.add(node.rest)
    return names


def _elements(target):
    if isinstance(target, (ast.Tuple, ast.List)):
        for element in target.elts:
            yield from _elements(element)
    elif isinstance(target, ast.Starred):
        yield from _elements(target.value)
    else:
        yield target


def _base(node):
    while isinstance(node, (ast.Subscript, ast.Attribute)):
        node = node.value
    return node.id if isinstance(node, ast.Name) else None


def _changed(stmts):
    """The names whose objects `stmts` change in place: an item or an attribute assigned or
    deleted, an augmented assignment, a method that changes a container."""
    names = set()
    for node in _nodes(stmts, comprehensions=False):
        targets = []
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, (ast.AugAssign, ast.AnnAssign, ast.For)):
            targets = [node.target]
        elif isinstance(node, ast.Delete):
            targets = node.targets
        elif isinstance(node, ast.withitem) and node.optional_vars is not None:
            targets = [node.optional_vars]
        for target in targets:
            for element in _elements(target):
                if isinstance(element, (ast.Subscript, ast.Attribute)):
                    names.add(_base(element))
                elif isinstance(element, ast.Name) and isinstance(node, ast.AugAssign):
                    names.add(element.id)
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in _MUTATORS
        ):
            names.add(_base(node.func.value))
    names.discard(None)
    return {name for name in names if not name.startswith(_PREFIX)}


def _parameters(arguments):
    names = [
        a.arg for a in arguments.posonlyargs + arguments.args + arguments.kwonlyargs
    ]
    names += [a.arg for a in (arguments.vararg, arguments.kwarg) if a is not None]
    return set(names)


def _generator(node):
    return any(
        isinstance(n, (ast.Yield, ast.YieldFrom, ast.Await)) for n in _nodes(node.body)
    )


def _branching(node):
    """Whether `node` has anything to convert."""
    for n in ast.walk(node):
        if isinstance(n, (ast.If, ast.IfExp, ast.BoolOp, ast.Assert, ast.Raise)):
            return True
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.Not):
            return True
        if isinstance(n, ast.Compare) and len(n.ops) > 1:
            return True
        if (
            isinstance(n, ast.Call)
            and isinstance(n.func, ast.Name)
            and n.func.id in _EXTREMA
        ):
            return True
    return False


# --------------------------------------------------------------------- building the tree


def _load(name):
    return ast.Name(id=name, ctx=ast.Load())


def _store(name, value):
    return ast.Assign(targets=[ast.Name(id=name, ctx=ast.Store())], value=value)


def _rt(attribute):
    return ast.Attribute(value=_load(_RT), attr=attribute, ctx=ast.Load())


def _call(function, arguments, keywords=()):
    return ast.Call(func=function, args=list(arguments), keywords=list(keywords))


def _tuple(elements):
    return ast.Tuple(elts=list(elements), ctx=ast.Load())


def _arguments(names=()):
    return ast.arguments(
        posonlyargs=[],
        args=[ast.arg(arg=name) for name in names],
        vararg=None,
        kwonlyargs=[],
        kw_defaults=[],
        kwarg=None,
        defaults=[],
    )


def _define(name, body, parameters=()):
    fields = {
        "name": name,
        "args": _arguments(parameters),
        "body": body,
        "decorator_list": [],
        "returns": None,
    }
    if "type_params" in ast.FunctionDef._fields:
        fields["type_params"] = []
    return ast.FunctionDef(**fields)


def _thunk(expression):
    for node in ast.walk(expression):
        if isinstance(node, ast.NamedExpr):
            raise Unconvertible(
                "an assignment expression (:=) in an operand of and, or, a conditional "
                "expression or a chained comparison, which runs in a function of its own",
                node,
            )
    return _point(ast.Lambda(args=_arguments(), body=expression), expression)


def _point(node, at):
    """Places `node` and the nodes in it that the conversion made (those without a place) at the
    start of `at`, the node of the source they come from: where a trace says they are. (A span
    over lines would put a call on its last line.)"""
    if getattr(at, "lineno", None) is None:
        return node
    if "lineno" in node._attributes and getattr(node, "lineno", None) is None:
        node.lineno = node.end_lineno = at.lineno
        node.col_offset = node.end_col_offset = at.col_offset
    for child in ast.iter_child_nodes(node):
        if "lineno" not in child._attributes or getattr(child, "lineno", None) is None:
            _point(child, at)
    return node


def _located(nodes, at):
    return [_point(node, at) for node in nodes]


class _Expressions(ast.NodeTransformer):
    """The expressions of a statement converted: a conditional expression, and, or, not, a
    chained comparison, min(), max(), any() and all()."""

    def visit_IfExp(self, node):
        self.generic_visit(node)
        call = _call(_rt("ifexp"), [node.test, _thunk(node.body), _thunk(node.orelse)])
        return _point(call, node)

    def visit_BoolOp(self, node):
        self.generic_visit(node)
        function = "and_" if isinstance(node.op, ast.And) else "or_"
        arguments = [node.values[0]] + [_thunk(value) for value in node.values[1:]]
        return _point(_call(_rt(function), arguments), node)

    def visit_UnaryOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.Not):
            return _point(_call(_rt("not_"), [node.operand]), node)
        return node

    def visit_Compare(self, node):
        self.generic_visit(node)
        if len(node.ops) == 1:
            return node
        arguments = [node.left]
        for op, comparator in zip(node.ops, node.comparators):
            arguments += [ast.Constant(type(op).__name__), _thunk(comparator)]
        return _point(_call(_rt("compare"), arguments), node)

    def visit_Call(self, node):
        self.generic_visit(node)
        if isinstance(node.func, ast.Name) and node.func.id in _EXTREMA:
            call = _call(
                _rt(f"call_{node.func.id}"), [node.func] + node.args, node.keywords
            )
            return _point(call, node)
        return node


# ---------------------------------------------------------------------- the conversion


class _Conversion:
    def __init__(self):
        self.count = 0

    def number(self):
        self.count += 1
        return self.count


class _Function:
    """The conversion of a function (or a lambda) of the model."""

    def __init__(self, conversion, node):
        self.conversion = conversion
        self.node = node
        self.expressions = _Expressions()
        self.parameters = _parameters(node.args)
        # the names a branch binds, which the function binds first: nonlocal in the branch
        self.unbound = set()
        # the loops around, innermost last: the names of their flags (broken, continued), or
        # None for a loop that breaks and continues as Python does
        self.loops = []
        self.globals = set()
        self.nonlocals = set()
        self.rewriting = False
        if isinstance(node, ast.Lambda):
            return
        for stmt in _statements(node.body):
            if isinstance(stmt, ast.Global):
                self.globals.update(stmt.names)
            elif isinstance(stmt, ast.Nonlocal):
                self.nonlocals.update(stmt.names)
        # a return in an if or a loop: every return sets the value and the flag that the
        # function has returned, and the rest of it runs where it has not
        self.rewriting = any(
            _has(ast.Return, block)
            for stmt in _statements(node.body)
            if isinstance(stmt, (ast.If,) + _LOOPS)
            for block in _blocks(stmt)
        )

    def expression(self, node):
        return self.expressions.visit(node) if node is not None else None

    def convert(self):
        node = self.node
        if isinstance(node, ast.Lambda):
            node.body = self.expression(node.body)
            return node
        body = node.body
        head = []
        if (
            body
            and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)
        ):
            head, body = body[:1], body[1:]
        leaves = _leaves(body)
        converted = self.block(body)
        prologue = []
        if self.rewriting:
            prologue += [
                _store(_RETURNED, ast.Constant(False)),
                _store(_RETVAL, _rt("UNDEFINED")),
            ]
        for name in sorted(self.unbound):
            prologue.append(_store(name, _call(_rt("unbound"), [ast.Constant(name)])))
        epilogue = []
        if self.rewriting:
            if leaves:
                value = _load(_RETVAL)
            else:
                value = _call(
                    _rt("result"),
                    [_load(_RETURNED), _load(_RETVAL), ast.Constant(node.name)],
                )
            epilogue.append(ast.Return(value))
        if body:
            _located(prologue, body[0])
            _located(epilogue, body[-1])
        node.body = head + prologue + converted + epilogue or [ast.Pass()]
        return node

    # ------------------------------------------------------------------------- blocks

    def block(self, stmts):
        """`stmts` converted: after a statement that returns, breaks or continues at some points,
        the rest runs where no point has."""
        converted = []
        for index, stmt in enumerate(stmts):
            # of the source: the conversion changes a statement in place
            leaves = _leaves([stmt])
            flags = self.flags(stmt)
            converted.extend(_located(self.statement(stmt), stmt))
            if leaves:
                break
            rest = stmts[index + 1 :]
            if flags and rest:
                converted.extend(self.guard(rest, flags))
                break
        return converted

    def flags(self, stmt):
        """The flags `stmt` may set: that the function returned, that its loop broke or went on
        with the next item."""
        flags = []
        if self.rewriting and _has(ast.Return, [stmt]):
            flags.append(_RETURNED)
        loop = self.loops[-1] if self.loops else None
        if loop is not None:
            broken, continued = loop
            if broken and _has(ast.Break, [stmt], loops=False):
                flags.append(broken)
            if continued and _has(ast.Continue, [stmt], loops=False):
                flags.append(continued)
        return flags

    def bound(self, stmts, loop=None):
        """The names a function of `stmts` binds: theirs, and the flags they set."""
        names = _bound(stmts)
        if self.rewriting and _has(ast.Return, stmts):
            names |= {_RETVAL, _RETURNED}
        if loop is None and self.loops:
            loop = self.loops[-1]
        if loop is not None:
            broken, continued = loop
            if broken and _has(ast.Break, stmts, loops=False):
                names.add(broken)
            if continued and _has(ast.Continue, stmts, loops=False):
                names.add(continued)
        return names

    def declarations(self, bound):
        self.unbound |= {
            name
            for name in bound
            if not name.startswith(_PREFIX)
            and name not in self.parameters | self.globals | self.nonlocals
        }
        declarations = []
        names = sorted(name for name in bound if name in self.globals)
        if names:
            declarations.append(ast.Global(names=names))
        names = sorted(name for name in bound if name not in self.globals)
        if names:
            declarations.append(ast.Nonlocal(names=names))
        return declarations

    def closure(self, name, stmts, bound, converted=None):
        if converted is None:
            converted = self.block(stmts)
        return _define(name, self.declarations(bound) + (converted or [ast.Pass()]))

    def accessors(self, number, bound, changed):
        """The functions that get the names a branch binds or changes, and set those it binds."""
        names = sorted(set(bound) | set(changed))
        get = _define(
            f"_sd__get{number}",
            [
                ast.Return(
                    ast.Dict(
                        keys=[ast.Constant(name) for name in names],
                        values=[_load(name) for name in names],
                    )
                )
            ],
        )
        values = "_sd__values"
        stores = [
            _store(
                name,
                ast.Subscript(
                    value=_load(values), slice=ast.Constant(name), ctx=ast.Load()
                ),
            )
            for name in sorted(bound)
        ]
        set_ = _define(
            f"_sd__set{number}",
            self.declarations(bound) + (stores or [ast.Pass()]),
            [values],
        )
        return get, set_

    def guard(self, stmts, flags, converted=None, bound=None):
        """`stmts` where none of `flags` holds."""
        number = self.conversion.number()
        if bound is None:
            bound = self.bound(stmts)
        changed = sorted(_changed(stmts))
        rest = self.closure(f"_sd__rest{number}", stmts, bound, converted)
        get, set_ = self.accessors(number, bound, changed)
        call = _call(
            _rt("unless"),
            [
                _tuple(_load(flag) for flag in flags),
                _load(rest.name),
                _load(get.name),
                _load(set_.name),
                _tuple(ast.Constant(name) for name in changed),
            ],
        )
        return _located([rest, get, set_, ast.Expr(call)], stmts[0])

    # --------------------------------------------------------------------- statements

    def statement(self, stmt):
        handler = getattr(self, f"statement_{type(stmt).__name__}", None)
        if handler is not None:
            return handler(stmt)
        if isinstance(stmt, (ast.AsyncFunctionDef, ast.ClassDef)):
            return [stmt]
        return [self.expressions.visit(stmt)]

    def statement_If(self, stmt):
        test = self.expression(stmt.test)
        number = self.conversion.number()
        stmts = stmt.body + stmt.orelse
        bound = self.bound(stmts)
        changed = sorted(_changed(stmts))
        then = self.closure(f"_sd__then{number}", stmt.body, bound)
        otherwise = self.closure(f"_sd__else{number}", stmt.orelse, bound)
        get, set_ = self.accessors(number, bound, changed)
        call = _call(
            _rt("if_stmt"),
            [
                test,
                _load(then.name),
                _load(otherwise.name),
                _load(get.name),
                _load(set_.name),
                _tuple(ast.Constant(name) for name in changed),
            ],
        )
        return [then, otherwise, get, set_, ast.Expr(call)]

    def statement_For(self, stmt):
        number = self.conversion.number()
        returns = self.rewriting and _has(ast.Return, stmt.body)
        breaks = _has(ast.Break, stmt.body, loops=False)
        continues = _has(ast.Continue, stmt.body, loops=False)
        # a break or a continue in a branch: flags, which every point of the loop has its own of
        flagged = returns or any(
            _has((ast.Break, ast.Continue), block, loops=False)
            for inner in _statements(stmt.body, loops=False)
            if isinstance(inner, ast.If)
            for block in _blocks(inner)
        )
        if isinstance(stmt, ast.While):
            stmt.test = self.expression(stmt.test)
        else:
            stmt.iter = self.expression(stmt.iter)
        if not flagged:
            self.loops.append(None)
            stmt.body = self.block(stmt.body) or [ast.Pass()]
            self.loops.pop()
            stmt.orelse = self.block(stmt.orelse)
            return [stmt]
        broken = f"_sd__broken{number}" if breaks else None
        continued = f"_sd__continued{number}" if continues else None
        flags = (broken, continued)
        bound = self.bound(stmt.body, loop=flags)
        self.loops.append(flags)
        body = self.block(stmt.body)
        self.loops.pop()
        # what carries over to the next items: a point that broke or returned stays out
        carried = [flag for flag in (broken, _RETURNED if returns else None) if flag]
        if carried:
            body = self.guard(stmt.body, carried, converted=body, bound=bound)
            if isinstance(stmt, ast.While):
                stmt.test = _call(
                    _rt("going"),
                    [_tuple(_load(flag) for flag in carried), _thunk(stmt.test)],
                )
            else:
                stmt.iter = _call(
                    _rt("loop"), [stmt.iter, _thunk(_tuple(_load(f) for f in carried))]
                )
        if continued:
            body = [_store(continued, ast.Constant(False))] + body
        stmt.body = body or [ast.Pass()]
        if stmt.orelse:
            # the else of a loop runs where it did not break
            if carried:
                stmt.orelse = self.guard(stmt.orelse, carried)
            else:
                stmt.orelse = self.block(stmt.orelse)
        before = [_store(broken, ast.Constant(False))] if broken else []
        return _located(before + [stmt], stmt)

    statement_While = statement_For

    def statement_Break(self, stmt):
        loop = self.loops[-1] if self.loops else None
        if loop is not None and loop[0]:
            return [_store(loop[0], ast.Constant(True))]
        return [stmt]

    def statement_Continue(self, stmt):
        loop = self.loops[-1] if self.loops else None
        if loop is not None and loop[1]:
            return [_store(loop[1], ast.Constant(True))]
        return [stmt]

    def statement_Return(self, stmt):
        value = self.expression(stmt.value)
        if not self.rewriting:
            stmt.value = value
            return [stmt]
        return [
            _store(_RETVAL, value if value is not None else ast.Constant(None)),
            _store(_RETURNED, ast.Constant(True)),
        ]

    def statement_Assert(self, stmt):
        source = ast.unparse(stmt.test)
        message = _thunk(self.expression(stmt.msg)) if stmt.msg is not None else None
        call = _call(
            _rt("assert_"),
            [
                self.expression(stmt.test),
                message if message is not None else ast.Constant(None),
                ast.Constant(source),
            ],
        )
        # as an assert, not run where Python runs with -O
        return [ast.If(test=_load("__debug__"), body=[ast.Expr(call)], orelse=[])]

    def statement_Raise(self, stmt):
        if stmt.exc is None:
            return [stmt]
        exception = "_sd__exception"
        raised = ast.Raise(exc=_load(exception), cause=self.expression(stmt.cause))
        return [
            _store(exception, self.expression(stmt.exc)),
            ast.If(
                test=_call(_rt("raising"), [_load(exception)]), body=[raised], orelse=[]
            ),
        ]

    def statement_Delete(self, stmt):
        converted = []
        others = []
        for target in stmt.targets:
            for element in _elements(target):
                if isinstance(element, ast.Name) and element.id not in self.globals:
                    if others:
                        converted.append(ast.Delete(targets=others))
                        others = []
                    # a name deleted has no value, as one not yet assigned
                    converted.append(
                        _store(
                            element.id,
                            _call(_rt("unbound"), [ast.Constant(element.id)]),
                        )
                    )
                else:
                    others.append(self.expression(element))
        if others:
            converted.append(ast.Delete(targets=others))
        return converted

    def statement_FunctionDef(self, stmt):
        stmt.decorator_list = [self.expression(d) for d in stmt.decorator_list]
        stmt.args.defaults = [self.expression(d) for d in stmt.args.defaults]
        stmt.args.kw_defaults = [self.expression(d) for d in stmt.args.kw_defaults]
        if _generator(stmt):
            return [stmt]
        return [_Function(self.conversion, stmt).convert()]

    def statement_With(self, stmt):
        for item in stmt.items:
            item.context_expr = self.expression(item.context_expr)
        stmt.body = self.block(stmt.body) or [ast.Pass()]
        return [stmt]

    def statement_Try(self, stmt):
        loop = self.loops[-1] if self.loops else None
        if (self.rewriting and _has(ast.Return, stmt.finalbody)) or (
            loop is not None
            and _has((ast.Break, ast.Continue), stmt.finalbody, loops=False)
        ):
            raise Unconvertible("a return, a break or a continue in a finally", stmt)
        stmt.body = self.block(stmt.body) or [ast.Pass()]
        for handler in stmt.handlers:
            handler.type = self.expression(handler.type)
            handler.body = self.block(handler.body) or [ast.Pass()]
        stmt.orelse = self.block(stmt.orelse)
        stmt.finalbody = self.block(stmt.finalbody)
        return [stmt]

    statement_TryStar = statement_Try

    def statement_Match(self, stmt):
        stmt.subject = self.expression(stmt.subject)
        for case in stmt.cases:
            case.guard = self.expression(case.guard)
            case.body = self.block(case.body) or [ast.Pass()]
        return [stmt]


# ------------------------------------------------------------------- source and code


def _source(function):
    """The syntax tree of the source of `function` at its lines, and its node in it; None where
    there is no source, or it is not one function of it."""
    code = function.__code__
    try:
        lines, first = inspect.getsourcelines(function)
    except (OSError, TypeError):
        return None
    if not lines:
        return None
    indented = lines[0][:1] in (" ", "\t")
    text = ("if 1:\n" if indented else "") + "".join(lines)
    try:
        tree = ast.parse(text, code.co_filename)
    except SyntaxError:
        return None
    ast.increment_lineno(tree, first - 2 if indented else first - 1)
    if code.co_name == "<lambda>":
        parameters = list(code.co_varnames[: code.co_argcount])
        found = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Lambda)
            and node.lineno == code.co_firstlineno
            and [a.arg for a in node.args.posonlyargs + node.args.args] == parameters
        ]
        return (tree, found[0]) if len(found) == 1 else None
    body = tree.body[0].body if indented else tree.body
    node = body[0] if body else None
    if not isinstance(node, ast.FunctionDef) or node.name != code.co_name:
        return None
    return tree, node


def _features(code):
    """The flags of the features of __future__ the code was compiled with."""
    flags = 0
    for name in __future__.all_feature_names:
        feature = getattr(__future__, name)
        mandatory = feature.getMandatoryRelease()
        if mandatory is not None and mandatory <= sys.version_info:
            continue
        if code.co_flags & feature.compiler_flag:
            flags |= feature.compiler_flag
    return flags


def _find(code, name, line):
    for constant in code.co_consts:
        if isinstance(constant, types.CodeType):
            if constant.co_name == name and constant.co_firstlineno == line:
                return constant
            found = _find(constant, name, line)
            if found is not None:
                return found
    return None


def _compile(node, code):
    """The code of `node`, a function or a lambda, compiled in a function that binds the free
    names of `code` -- and the runtime of the branches -- as the closure of the new function
    gives them; at the lines and with the features of `code`."""
    body = [_store(_RT, ast.Constant(None))]
    body += [_store(name, ast.Constant(None)) for name in code.co_freevars]
    if isinstance(node, ast.Lambda):
        body.append(ast.Return(node))
    else:
        body += [node, ast.Return(_load(node.name))]
    outer = _point(_define("_sd__outer", body), node)
    module = ast.Module(body=[outer], type_ignores=[])
    ast.fix_missing_locations(module)
    compiled = compile(
        module, code.co_filename, "exec", flags=_features(code), dont_inherit=True
    )
    return _find(compiled, code.co_name, code.co_firstlineno)


def _same_parameters(node, code):
    """Whether the function `node` of the source takes the parameters of `code`: a source that
    has changed since may not be called as the function is."""
    arguments = node.args
    names = [
        a.arg for a in arguments.posonlyargs + arguments.args + arguments.kwonlyargs
    ]
    for extra, flag in (
        (arguments.vararg, inspect.CO_VARARGS),
        (arguments.kwarg, inspect.CO_VARKEYWORDS),
    ):
        if (extra is not None) != bool(code.co_flags & flag):
            return False
        if extra is not None:
            names.append(extra.arg)
    count = code.co_argcount + code.co_kwonlyargcount
    count += bool(code.co_flags & inspect.CO_VARARGS) + bool(
        code.co_flags & inspect.CO_VARKEYWORDS
    )
    return (
        len(arguments.posonlyargs) == code.co_posonlyargcount
        and len(arguments.kwonlyargs) == code.co_kwonlyargcount
        and names == list(code.co_varnames[:count])
    )


def _warn(function, reason, line=None):
    code = function.__code__
    warnings.warn_explicit(
        f"{function.__qualname__} runs as it is written, without its branches on values of the "
        f"point converted: {reason}",
        UserWarning,
        code.co_filename,
        line if line is not None else code.co_firstlineno,
    )


def convert(function):
    """`function` with its branches converted: a new function of its globals, its closure and
    its defaults. None where it has nothing to convert, is not a function of Python with its
    source (a builtin, a generator, a wrapper), or the conversion cannot follow it (with a
    warning)."""
    if not isinstance(function, types.FunctionType) or hasattr(function, "__wrapped__"):
        return None
    code = function.__code__
    if code.co_flags & (
        inspect.CO_GENERATOR | inspect.CO_COROUTINE | inspect.CO_ASYNC_GENERATOR
    ):
        return None
    if "__class__" in code.co_freevars:
        return None
    found = _source(function)
    if found is None:
        return None
    _, node = found
    if not _branching(node):
        return None
    # the source as it is now, which a core file records; it takes the parameters of the function
    if not _same_parameters(node, code):
        _warn(function, "its source has changed since it was imported")
        return None
    try:
        converted = _Function(_Conversion(), node).convert()
        new_code = _compile(converted, code)
    except Unconvertible as error:
        _warn(function, str(error), getattr(error.node, "lineno", None))
        return None
    if new_code is None:
        return None
    cells = []
    for name in new_code.co_freevars:
        if name == _RT:
            cells.append(types.CellType(_branches))
        elif name in code.co_freevars:
            cells.append(function.__closure__[code.co_freevars.index(name)])
        else:
            return None
    new = types.FunctionType(
        new_code,
        function.__globals__,
        function.__name__,
        function.__defaults__,
        tuple(cells) if cells else None,
    )
    new.__kwdefaults__ = function.__kwdefaults__
    new.__qualname__ = function.__qualname__
    new.__module__ = function.__module__
    new.__doc__ = function.__doc__
    new.__dict__.update(function.__dict__)
    new.__wrapped__ = function
    return new
