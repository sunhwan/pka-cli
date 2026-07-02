#!/usr/bin/env bash
# Clones MolGpKa at the pinned SHA used across reactable.net and
# wash-benchmark, so all three stay reproducible against the same model.
# Also finishes installing torch-scatter/torch-sparse, which environment.yml
# deliberately leaves out (see the comment there) because they need staged,
# non-isolated installation.
set -euo pipefail

MOLGPKA_SHA="4dc8352"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="$REPO_ROOT/vendor/MolGpKa"

if [ -d "$DEST" ]; then
    echo "MolGpKa already present at $DEST"
else
    mkdir -p "$REPO_ROOT/vendor"
    git clone https://github.com/Xundrug/MolGpKa.git "$DEST"
    git -C "$DEST" checkout "$MOLGPKA_SHA"
    echo "MolGpKa checked out at $MOLGPKA_SHA -> $DEST"
fi

if python -c "import torch_scatter, torch_sparse" >/dev/null 2>&1; then
    echo "torch-scatter/torch-sparse already installed"
    exit 0
fi

echo "Installing torch-scatter/torch-sparse (no prebuilt wheels for this platform; building from source)..."

# torch-scatter/torch-sparse `import torch` inside their own setup.py to
# build CPU extensions against libtorch, so they must be installed with
# --no-build-isolation (after torch is already present) or the build can't
# see torch at all.
#
# Separately, PyTorch 2.2.2's bundled headers (c10/util/strong_type.h)
# specialize std::is_arithmetic, which is ill-formed C++; newer clang
# (Xcode 16+, and recent conda-forge toolchains) rejects this as a hard
# error (-Winvalid-specialization) instead of a warning, so any C++
# extension built against these torch headers fails to compile. This is a
# bug in PyTorch's headers, not in torch-scatter/torch-sparse, so we just
# downgrade the diagnostic back to non-fatal.
CXXFLAGS="${CXXFLAGS:-} -Wno-invalid-specialization" \
CFLAGS="${CFLAGS:-} -Wno-invalid-specialization" \
    pip install --no-build-isolation \
    torch-scatter==2.1.2 torch-sparse==0.6.18 \
    --extra-index-url https://data.pyg.org/whl/torch-2.2.0+cpu.html

echo "torch-scatter/torch-sparse installed"
