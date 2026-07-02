# pka-cli

Standalone pKa site prediction + protomer enumeration, extracted from
[reactable.net](../reactable.net)'s Modal pKa backend. Wraps
[MolGpKa](https://github.com/Xundrug/MolGpKa) (a CPU GCN predicting
per-ionizable-site pKa) with an independent-site microstate enumeration
model, as a local CLI/library instead of a hosted endpoint.

Includes two fixes not present in upstream MolGpKa or its naive usage:
1. Acid sites are indexed by MolGpKa as the acidic *hydrogen* atom, not the
   heavy atom bearing it — mapped to the heavy-atom neighbor so results
   aren't silently dropped or mismatched against base sites.
2. Symmetric groups (e.g. primary `-NH2`) produce duplicate acid sites
   (one per equivalent H) — deduped by `(kind, atom_idx)`.

## Setup

```bash
conda env create -f environment.yml
conda activate pka-cli
./scripts/setup_molgpka.sh   # clones MolGpKa at the pinned SHA into vendor/
pip install -e .
```

## Usage

```bash
pka-cli "CC(=O)Nc1ccc(O)cc1" --ph 7.4 --max-protomers 8
```

```json
{
  "CC(=O)Nc1ccc(O)cc1": {
    "neutral_smiles": "CC(=O)Nc1ccc(O)cc1",
    "sites": [
      {"kind": "acid", "atom_idx": 8, "pka": 10.1},
      {"kind": "acid", "atom_idx": 3, "pka": 13.84},
      {"kind": "base", "atom_idx": 3, "pka": 4.6}
    ],
    "microspecies": [
      {"smiles": "CC(=O)Nc1ccc(O)cc1", "charge": 0, "population": 0.9964},
      {"smiles": "CC(=O)Nc1ccc([O-])cc1", "charge": -1, "population": 0.002},
      {"smiles": "CC(=O)[NH2+]c1ccc(O)cc1", "charge": 1, "population": 0.0016}
    ],
    "dominant_pka": 10.1,
    "ph": 7.4
  }
}
```

Read from a file (one SMILES per line) with `--input molecules.smi`, or use
`pka_cli.core.predict_batch(...)` directly as a library.

## Tests

```bash
pytest
```

Tests load the real MolGpKa models and are skipped automatically if
`vendor/MolGpKa` hasn't been set up yet.
