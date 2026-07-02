#!/usr/bin/env bash
# Clones MolGpKa at the pinned SHA used across reactable.net and
# wash-benchmark, so all three stay reproducible against the same model.
set -euo pipefail

MOLGPKA_SHA="4dc8352"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="$REPO_ROOT/vendor/MolGpKa"

if [ -d "$DEST" ]; then
    echo "MolGpKa already present at $DEST"
    exit 0
fi

mkdir -p "$REPO_ROOT/vendor"
git clone https://github.com/Xundrug/MolGpKa.git "$DEST"
git -C "$DEST" checkout "$MOLGPKA_SHA"
echo "MolGpKa checked out at $MOLGPKA_SHA -> $DEST"
