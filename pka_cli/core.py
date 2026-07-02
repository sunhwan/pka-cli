"""Standalone MolGpKa pKa prediction + protomer enumeration.

Ported from reactable.net's backend/modal_app.py, with two fixes applied
upstream there and carried over here:

1. MolGpKa indexes *acid* ionization sites by the acidic hydrogen atom
   itself, not the heavy atom bearing it (base sites are already
   heavy-atom indexed). Naively filtering "keep only heavy-atom indices"
   silently drops every acid site. We map H -> its heavy-atom neighbor.
2. Symmetric groups (e.g. primary -NH2) have multiple equivalent
   ionizable H's, so MolGpKa's site detector returns one aid per H. After
   the H -> heavy-atom mapping these collide on the same atom_idx; we
   dedupe by (kind, atom_idx), averaging pKa across the collision.
"""
import itertools
import os.path as osp
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from . import _env


@dataclass
class Site:
    kind: str  # "acid" or "base"
    atom_idx: int  # heavy-atom index into the neutral (no explicit-H) mol
    pka: float


@dataclass
class Microspecies:
    smiles: str
    charge: int
    population: float


@dataclass
class PkaResult:
    smiles: str
    neutral_smiles: str
    sites: List[Site]
    microspecies: List[Microspecies]
    dominant_pka: Optional[float]
    ph: float


class PkaPredictor:
    """Loads MolGpKa's two GCN models once and reuses them across calls."""

    def __init__(self, molgpka_home: Optional[str] = None):
        if molgpka_home is not None:
            import os
            os.environ["MOLGPKA_HOME"] = molgpka_home
        self._home, self._src = _env.ensure_importable()
        self._acid_model, self._base_model = self._load_models()

    def _load_models(self):
        import torch
        from utils.net import GCNNet

        def _load(path):
            m = GCNNet().to("cpu")
            m.load_state_dict(torch.load(path, map_location="cpu"))
            m.eval()
            return m

        models_dir = osp.join(self._src, "..", "models")
        acid = _load(osp.join(models_dir, "weight_acid.pth"))
        base = _load(osp.join(models_dir, "weight_base.pth"))
        return acid, base

    def _predict_sites(self, smiles: str) -> Tuple[Optional[List[Site]], Optional[str]]:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.MolStandardize import rdMolStandardize
        import torch
        from utils.descriptor import mol2vec
        from utils.ionization_group import get_ionization_aid

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        neutral_smiles = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(neutral_smiles)
        if mol is None:
            return None, None
        # AddHs appends H atoms after heavy atoms, so heavy-atom indices are
        # stable between the H-added mol (used for the GCN) and the neutral
        # SMILES mol.
        molh = AllChem.AddHs(mol)

        def _pka(aid, model):
            data = mol2vec(molh, aid)
            with torch.no_grad():
                out = model(data.to("cpu"))
            return float(out.cpu().numpy()[0][0])

        def _is_heavy(aid):
            return 0 <= aid < molh.GetNumAtoms() and molh.GetAtomWithIdx(aid).GetAtomicNum() != 1

        def _heavy_parent(aid):
            if not (0 <= aid < molh.GetNumAtoms()):
                return None
            atom = molh.GetAtomWithIdx(aid)
            if atom.GetAtomicNum() != 1:
                return aid
            nbrs = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1]
            return nbrs[0] if nbrs else None

        by_key: Dict[Tuple[str, int], List[float]] = {}
        for aid in get_ionization_aid(molh, acid_or_base="acid"):
            heavy_idx = _heavy_parent(aid)
            if heavy_idx is not None:
                by_key.setdefault(("acid", int(heavy_idx)), []).append(_pka(aid, self._acid_model))
        for aid in get_ionization_aid(molh, acid_or_base="base"):
            if _is_heavy(aid):
                by_key.setdefault(("base", int(aid)), []).append(_pka(aid, self._base_model))

        sites = [Site(kind=kind, atom_idx=atom_idx, pka=sum(pkas) / len(pkas))
                 for (kind, atom_idx), pkas in by_key.items()]
        return sites, neutral_smiles

    @staticmethod
    def _build_microstate_smiles(neutral_smiles, site_states):
        """site_states: list of (atom_idx, kind, protonated). Returns (smiles, charge)."""
        from rdkit import Chem

        base = Chem.MolFromSmiles(neutral_smiles)
        if base is None:
            return None, 0
        molh = Chem.AddHs(base)
        rw = Chem.RWMol(molh)
        n = rw.GetNumAtoms()
        charge = 0
        to_remove = []
        for atom_idx, kind, protonated in site_states:
            if atom_idx < 0 or atom_idx >= n:
                return None, charge
            atom = rw.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 1:  # never edit an H "site"
                continue
            if kind == "base":
                if protonated:  # add a proton -> +1
                    atom.SetFormalCharge(1)
                    h = rw.AddAtom(Chem.Atom(1))
                    rw.AddBond(atom_idx, h, Chem.BondType.SINGLE)
                    charge += 1
            else:  # acid
                if not protonated:  # remove a proton -> -1
                    h_nbr = next((nb.GetIdx() for nb in atom.GetNeighbors()
                                  if nb.GetAtomicNum() == 1), None)
                    if h_nbr is None:
                        return None, charge
                    atom.SetFormalCharge(-1)
                    to_remove.append(h_nbr)
                    charge -= 1
        for idx in sorted(set(to_remove), reverse=True):
            rw.RemoveAtom(idx)
        try:
            m = rw.GetMol()
            Chem.SanitizeMol(m)
            m = Chem.RemoveHs(m)
            return Chem.MolToSmiles(m), charge
        except Exception:
            return None, charge

    def _enumerate_microspecies(self, neutral_smiles, sites: List[Site], ph, max_protomers):
        # Cap site count to keep 2^n manageable; keep sites nearest pH.
        usable = sorted(sites, key=lambda s: abs(s.pka - ph))[:8]
        if not usable:
            return [Microspecies(smiles=neutral_smiles, charge=0, population=1.0)]

        def p_prot(pka):
            return 1.0 / (1.0 + 10.0 ** (ph - pka))

        fracs = [p_prot(s.pka) for s in usable]
        out: Dict[str, Microspecies] = {}
        for combo in itertools.product([True, False], repeat=len(usable)):
            pop = 1.0
            states = []
            for s, frac, protonated in zip(usable, fracs, combo):
                pop *= frac if protonated else (1.0 - frac)
                states.append((s.atom_idx, s.kind, protonated))
            smi, charge = self._build_microstate_smiles(neutral_smiles, states)
            if smi is None:
                continue
            if smi in out:
                out[smi].population += pop
            else:
                out[smi] = Microspecies(smiles=smi, charge=charge, population=pop)
        ranked = sorted(out.values(), key=lambda m: m.population, reverse=True)
        return ranked[:max_protomers]

    def predict(self, smiles: str, ph: float = 7.4, max_protomers: int = 5) -> Optional[PkaResult]:
        sites, neutral = self._predict_sites(smiles)
        if sites is None:
            return None
        micro = self._enumerate_microspecies(neutral, sites, ph, max_protomers)
        dominant = min(sites, key=lambda s: abs(s.pka - ph)).pka if sites else None
        return PkaResult(
            smiles=smiles,
            neutral_smiles=neutral,
            sites=sites,
            microspecies=micro,
            dominant_pka=dominant,
            ph=ph,
        )


def predict_batch(
    smiles_list: List[str],
    ph: float = 7.4,
    max_protomers: int = 5,
    molgpka_home: Optional[str] = None,
) -> Dict[str, Optional[PkaResult]]:
    """Convenience wrapper: loads the models once, predicts for every SMILES."""
    predictor = PkaPredictor(molgpka_home=molgpka_home)
    results = {}
    for smi in smiles_list:
        if not smi or smi in results:
            continue
        results[smi] = predictor.predict(smi, ph=ph, max_protomers=max_protomers)
    return results
