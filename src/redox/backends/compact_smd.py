"""Reuse solvent-build temporaries without changing SMD surface or equations."""

import ast
from contextlib import contextmanager
import inspect
import textwrap


def compact_build(function):
    """Transform only the installed dense K = S - factor*DAS assignment."""
    tree = ast.parse(textwrap.dedent(inspect.getsource(function)))
    changed = 0
    body = []
    for statement in tree.body[0].body:
        if (
            isinstance(statement, ast.Assign)
            and len(statement.targets) == 1
            and isinstance(statement.targets[0], ast.Name)
            and statement.targets[0].id == "K"
        ):
            expression = statement.value
            if isinstance(expression, ast.Constant) and expression.value is None:
                body.append(statement)
                body.extend(ast.parse("DAS = None").body)
                continue
            expected = ast.parse("S - f_epsilon/(2.0*np.pi) * DAS", mode="eval").body
            if ast.dump(expression) != ast.dump(expected):
                raise ValueError(
                    "Installed SMD matrix equation differs from the reviewed source"
                )
            body.extend(
                ast.parse(
                    "DA = None\nDAS *= -f_epsilon/(2.0*np.pi)\nDAS += S\nK = DAS"
                ).body
            )
            changed += 1
        else:
            body.append(statement)
    if changed != 1:
        raise ValueError("Expected exactly one reviewed SMD matrix assignment")
    tree.body[0].body = body
    ast.fix_missing_locations(tree)
    namespace = {}
    exec(compile(tree, "<compact_smd_build>", "exec"), function.__globals__, namespace)
    return namespace[function.__name__], ast.unparse(tree)


@contextmanager
def compact_solvent_build():
    from gpu4pyscf.solvent import smd

    original = smd.SMD.build
    replacement, source = compact_build(original)
    smd.SMD.build = replacement
    try:
        yield source
    finally:
        smd.SMD.build = original
