#!/usr/bin/env bash
# update_lock.sh - regenerate conda-lock-linux-64.yaml from conda.yaml
#
# input:   conda.yaml at the repo root
# output:  conda-lock-linux-64.yaml, overwritten in place
# usage:   scripts/update_lock.sh
# env:     conda-lock on PATH
# notes:   singularity/qutrna2.def installs the lockfile, and a conda.yaml
#          edit reaches the container only once both are committed. the name
#          avoids conda-lock.yaml, which conda routes to its own parser
set -euo pipefail

cd "$(dirname "$0")/.."
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

conda-lock lock -f conda.yaml -p linux-64 --lockfile "$tmp/lock.yml"
conda-lock render -k env --filename-template "$tmp/rendered-{platform}" "$tmp/lock.yml"

# sprinx comes from git and cannot be hashed. pip needs hashes on every
# requirement or on none.
sed 's/ --hash=sha256:[0-9a-f]*//' "$tmp/rendered-linux-64.yml" > conda-lock-linux-64.yaml

echo "wrote conda-lock-linux-64.yaml"
