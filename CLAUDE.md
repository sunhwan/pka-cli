# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`pka-cli` is a standalone pKa site prediction + protomer enumeration tool, extracted
from reactable.net's Modal pKa backend. It wraps
[MolGpKa](https://github.com/Xundrug/MolGpKa) (a CPU GCN that predicts per-ionizable-site
pKa) with an independent-site microstate enumeration model, exposed as a local
CLI/library instead of a hosted endpoint.

MolGpKa itself is *not* vendored in git — it's a pinned external checkout (see
Architecture below), so it won't show up in `git log`/`git blame` for this repo.

## Setup

```bash
conda env create -f environment.yml
conda activate pka-cli
./scripts/setup_molgpka.sh   # clones MolGpKa + finishes the torch-geometric extension install
pip install -e .
```

`environment.yml` intentionally does *not* list `torch-scatter`/`torch-sparse`: they
`import torch` inside their own `setup.py` to build CPU extensions, and pip's build
isolation hides an already-installed `torch` from that build step, so a single combined
`pip install` fails. `scripts/setup_molgpka.sh` installs them afterward, separately, with
`--no-build-isolation`, and also passes `-Wno-invalid-specialization` — PyTorch 2.2.2's
bundled headers do something newer clang (Xcode 16+, recent conda-forge toolchains)
treats as a hard error rather than a warning, which otherwise breaks the build. The
script is idempotent (no-ops the clone if `vendor/MolGpKa` already exists, and skips the
pip step if `torch_scatter`/`torch_sparse` already import) and pins MolGpKa to a specific
commit SHA so results stay reproducible against the same model weights across this repo,
reactable.net, and wash-benchmark.

## Commands

```bash
# Run the CLI
pka-cli "CC(=O)Nc1ccc(O)cc1" --ph 7.4 --max-protomers 8
pka-cli --input molecules.smi   # one SMILES per line

# Run all tests
pytest

# Run a single test
pytest tests/test_core.py::test_acetaminophen_dominant_site_is_the_phenol
```

Tests load the real MolGpKa GCN models via the session-scoped `predictor` fixture
(`tests/conftest.py`) and are auto-skipped if `vendor/MolGpKa` hasn't been set up.

## Architecture

Three-module core, deliberately kept small:

- **`pka_cli/_env.py`** — locates the MolGpKa checkout (`MOLGPKA_HOME` env var, or
  `vendor/MolGpKa` by default) and makes it importable. This is fiddly because MolGpKa
  is a foreign codebase, not a proper package: `ensure_importable()` inserts its `src/`
  onto `sys.path` *and* `chdir`s into it, because MolGpKa's ionization-site detector
  reads `utils/smarts_pattern.tsv` relative to cwd. It also stubs out
  `scipy.sparse.linalg._propack` on import failure — a workaround for a broken PROPACK
  Fortran extension on some macOS builds; MolGpKa's GCN inference never touches scipy's
  SVD solvers, so the stub is safe.

- **`pka_cli/core.py`** — the actual prediction logic, in `PkaPredictor`:
  1. Loads MolGpKa's two GCN weights (acid/base) once per instance.
  2. `_predict_sites`: standardizes/neutralizes the input SMILES with RDKit, calls into
     MolGpKa's `get_ionization_aid` + `mol2vec` to get per-site pKa. Two fixes are
     applied here that are *not* in upstream MolGpKa or naive usage of it — see the
     module docstring and README for the exact mechanics (H-atom-vs-heavy-atom acid site
     indexing, and dedup of symmetric-group duplicate sites). These fixes are the reason
     this repo exists rather than calling MolGpKa directly; don't regress them.
  3. `_enumerate_microspecies`: treats ionizable sites as independent (Henderson-
     Hasselbalch per site, no cross-site coupling), enumerates all 2^n protonation
     combinations (capped at the 8 sites nearest the target pH to keep this tractable),
     builds each microstate's SMILES via RDKit atom edits, and aggregates duplicate
     SMILES that different protonation combos collapse to.
  4. `predict_batch` is the library entry point: loads models once, predicts over a list
     of SMILES.

- **`pka_cli/cli.py`** — thin argparse wrapper over `predict_batch`, dumps JSON to
  stdout. SMILES come from positional args and/or `--input FILE` (one per line, `#`
  comments ignored).

### Data flow for a single SMILES

`SMILES -> RDKit neutralize -> AddHs (stable heavy-atom indices) -> MolGpKa site
detection (acid + base GCNs) -> Site list (deduped, heavy-atom indexed) -> per-site
protonation fractions at target pH -> enumerate microstates -> rank by population ->
PkaResult`

Heavy-atom indices are load-bearing across this pipeline: `AddHs` appends H atoms after
existing heavy atoms, so heavy-atom indices are guaranteed stable between the
H-added mol used for GCN inference and the neutral-SMILES mol used for reporting/microstate
building. Don't reorder atoms between these steps.
