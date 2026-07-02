"""Command-line interface: predict pKa sites and enumerate protomers."""
import argparse
import dataclasses
import json
import sys

from .core import predict_batch


def _read_smiles(args) -> list:
    smiles = list(args.smiles)
    if args.input:
        with open(args.input) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    smiles.append(line.split()[0])
    if not smiles:
        raise SystemExit("no SMILES given (pass positional args or --input FILE)")
    return smiles


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="pka-cli",
        description="Predict per-site pKa and enumerate protomers with MolGpKa.",
    )
    parser.add_argument("smiles", nargs="*", help="SMILES string(s) to predict")
    parser.add_argument("--input", "-i", help="file with one SMILES per line")
    parser.add_argument("--ph", type=float, default=7.4, help="target pH (default: 7.4)")
    parser.add_argument("--max-protomers", type=int, default=5, help="max microspecies to report")
    parser.add_argument("--molgpka-home", help="path to a MolGpKa checkout (overrides MOLGPKA_HOME)")
    parser.add_argument("--indent", type=int, default=2, help="JSON indent (0 for compact)")
    args = parser.parse_args(argv)

    smiles_list = _read_smiles(args)
    results = predict_batch(
        smiles_list,
        ph=args.ph,
        max_protomers=args.max_protomers,
        molgpka_home=args.molgpka_home,
    )

    out = {}
    for smi, result in results.items():
        if result is None:
            out[smi] = None
            continue
        out[smi] = {
            "neutral_smiles": result.neutral_smiles,
            "sites": [dataclasses.asdict(s) for s in result.sites],
            "microspecies": [dataclasses.asdict(m) for m in result.microspecies],
            "dominant_pka": result.dominant_pka,
            "ph": result.ph,
        }

    json.dump(out, sys.stdout, indent=args.indent or None)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
