"""test_rules.py - each subcommand named by a rule exists in its script

input:   workflow/rules/*.smk, workflow/Snakefile, workflow/scripts/*.py
output:  pytest assertions
usage:   pytest tests/test_rules.py
env:     none
notes:   the scripts are read with ast. The check therefore runs in any
         environment that carries pytest. A dry run builds the DAG but does
         not execute a rule body, and a wrong subcommand survives it
"""
import ast
import os
import re

import pytest

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RULES = os.path.join(REPO, "workflow", "rules")
SCRIPTS = os.path.join(REPO, "workflow", "scripts")

# the token that follows the script path in a rule
INVOCATION = re.compile(r"scripts/(\w+)\.py\s+([a-z][a-z0-9-]*)")


def rule_sources():
    yield os.path.join(REPO, "workflow", "Snakefile")
    for name in sorted(os.listdir(RULES)):
        if name.endswith(".smk"):
            yield os.path.join(RULES, name)


def command_name(decorator, function):
    """The name passed to click, or the function name with dashes."""
    if isinstance(decorator, ast.Call):
        for arg in decorator.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                return arg.value
        for kw in decorator.keywords:
            if kw.arg == "name" and isinstance(kw.value, ast.Constant):
                return kw.value.value
    return function.replace("_", "-")


def click_commands(script):
    """The subcommands defined by the script's click group.

    Returns None for a script without one.
    """
    with open(os.path.join(SCRIPTS, f"{script}.py")) as fh:
        tree = ast.parse(fh.read())
    group, names = False, set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        for decorator in node.decorator_list:
            text = ast.unparse(decorator)
            if text.startswith("click.group"):
                group = True
            elif text.startswith("cli.command"):
                names.add(command_name(decorator, node.name))
    return names if group else None


def invocations():
    """The (script, subcommand) pairs named by the workflow, click groups only."""
    found = set()
    for path in rule_sources():
        with open(path) as fh:
            # a rule may break the line before the subcommand
            text = fh.read().replace("\\\n", " ")
        found.update(INVOCATION.findall(text))
    return sorted(p for p in found if click_commands(p[0]) is not None)


PAIRS = invocations()


def test_the_scan_reaches_every_group_script():
    """An empty scan would make the checks below pass without testing anything."""
    groups = {name[:-3] for name in os.listdir(SCRIPTS)
              if name.endswith(".py") and click_commands(name[:-3]) is not None}
    assert {script for script, _ in PAIRS} == groups


@pytest.mark.parametrize("script,subcommand", PAIRS)
def test_rule_subcommand_exists(script, subcommand):
    defined = click_commands(script)
    assert subcommand in defined, \
        f"a rule runs {script}.py {subcommand}; it defines {sorted(defined)}"
