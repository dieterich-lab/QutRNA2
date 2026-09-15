"""test_packaging.py - static checks on the files installed by the container

input:   singularity/qutrna2.def, conda.yaml, the lockfile, the config schema
output:  pytest assertions
usage:   pytest tests/test_packaging.py
env:     pyyaml. QUTRNA2_REQUIRE_TOOLS=1 makes the binaries mandatory
notes:   apptainer does not offer a parse-only mode, and a hosted runner
         lacks the disk to build. The recipe is read as text
"""
import os
import re
import shutil
import subprocess

import pytest
import yaml

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

# conda routes these two names to its own lockfile parser
RESERVED = ("conda-lock.yml", "conda-lock.yaml")
# name=version only. A range does not pin a version that can be compared
PIN = re.compile(r"^([A-Za-z0-9_.-]+)=([0-9][^=]*)$")
INSTALL = re.compile(r"conda env create\s+-f\s+(\S+)")


def load(*parts):
    with open(os.path.join(REPO, *parts)) as fh:
        return yaml.safe_load(fh)


def recipe_section(name):
    """Non-blank stripped lines of one %section of qutrna2.def."""
    with open(os.path.join(REPO, "singularity", "qutrna2.def")) as fh:
        lines = fh.read().splitlines()
    out, inside = [], False
    for line in lines:
        if line.startswith("%"):
            inside = line.strip() == name
        elif inside and line.strip():
            out.append(line.strip())
    return out


def test_recipe_files_exist_and_feed_post():
    """A %files path resolves against the build directory."""
    build = os.path.join(REPO, "singularity")
    copied = set()
    for entry in recipe_section("%files"):
        source = entry.split()[0]
        assert os.path.exists(os.path.join(build, source)), f"%files names {source}"
        copied.add(os.path.basename(source))

    match = INSTALL.search(" ".join(recipe_section("%post")))
    assert match, "%post does not run 'conda env create -f'"
    named = os.path.basename(match.group(1))
    assert named in copied, f"%post installs {named}; %files copies {sorted(copied)}"
    assert named not in RESERVED, f"{named} reaches conda's lockfile parser"


def test_every_pin_is_locked_at_that_version():
    lock = {}
    for dep in load("conda-lock-linux-64.yaml")["dependencies"]:
        if isinstance(dep, str):
            name, _, rest = dep.partition("=")
            lock[name.lower()] = rest.partition("=")[0]

    for dep in load("conda.yaml")["dependencies"]:
        match = PIN.match(dep) if isinstance(dep, str) else None
        if not match:
            continue
        # conda lowercases names: JACUSA2 is locked as jacusa2
        name, version = match.group(1).lower(), match.group(2)
        assert name in lock, f"conda.yaml pins {name} and the lockfile omits it"
        # a pin may be shorter than the locked version: 9.8 covers 9.8.2
        depth = version.count(".") + 1
        assert lock[name].split(".")[:depth] == version.split("."), \
            (f"conda.yaml pins {name}={version}, lockfile carries {lock[name]}. "
             "Run scripts/update_lock.sh")


def test_sprinx_is_pinned_to_an_immutable_ref():
    """A floating ref resolves at install time and two builds then differ."""
    deps = load("conda.yaml")["dependencies"]
    pip = next(d["pip"] for d in deps if isinstance(d, dict) and "pip" in d)
    ref = next(p for p in pip if "sprinx" in p).rpartition("@")[2]
    assert ref not in ("main", "master", "HEAD"), f"sprinx tracks {ref}"


# read from the schema, which is overridden by a config
JACUSA2 = load("workflow", "schemas", "config.yaml")[
    "properties"]["jacusa2"]["properties"]["jar"]["default"]

# gpu-tRNA-mapper is left out; probing it needs a GPU
TOOLS = [(JACUSA2, "-h"), ("cmalign", "-h"), ("samtools", "--version")]


@pytest.mark.parametrize("binary,arg", TOOLS)
def test_tool_answers(binary, arg):
    """A missing binary surfaces as a bare 'command not found' inside a rule."""
    found = shutil.which(binary)
    if found is None:
        if os.environ.get("QUTRNA2_REQUIRE_TOOLS") == "1":
            pytest.fail(f"{binary} is absent from PATH")
        pytest.skip(f"{binary} is absent from PATH")
    proc = subprocess.run([found, arg], capture_output=True, text=True, check=False)
    assert proc.stdout or proc.stderr, f"{binary} did not answer"
