"""Validation against literature pKa values, mirroring the manual checks
run against the reactable.net production endpoint before/after the two
site-detection fixes (see core.py's module docstring).
"""
import pytest

ACETAMINOPHEN = "CC(=O)Nc1ccc(O)cc1"
PHENOL = "c1ccccc1O"
TERT_BUTYLPHENOL = "CC(C)(C)c1ccc(O)cc1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
BENZOIC_ACID = "c1ccccc1C(=O)O"
ANILINE = "c1ccccc1N"
PYRIDINE = "c1ccncc1"
IMIDAZOLE = "c1cnc[nH]1"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
SULFANILAMIDE = "Nc1ccc(cc1)S(=O)(=O)N"


def test_acetaminophen_dominant_site_is_the_phenol(predictor):
    """Regression test for the H-vs-heavy-atom acid site bug: before the
    fix, every acid site was silently dropped and the dominant site fell
    back to the amide N (pKa ~4.6), missing the phenol entirely."""
    result = predictor.predict(ACETAMINOPHEN, ph=7.4, max_protomers=8)
    assert result is not None
    acid_sites = [s for s in result.sites if s.kind == "acid"]
    assert acid_sites, "expected at least one acid site (the phenol)"
    assert result.dominant_pka == pytest.approx(10.1, abs=1.0)

    neutral_pop = next(
        (m.population for m in result.microspecies if m.charge == 0), None
    )
    assert neutral_pop is not None and neutral_pop > 0.9


@pytest.mark.parametrize(
    "smiles,expected_pka",
    [
        (PHENOL, 9.95),
        (TERT_BUTYLPHENOL, 10.23),
        (ASPIRIN, 3.5),
        (IBUPROFEN, 4.6),
        (BENZOIC_ACID, 4.2),
    ],
)
def test_acid_pka_matches_literature(predictor, smiles, expected_pka):
    result = predictor.predict(smiles, ph=7.4, max_protomers=5)
    assert result is not None
    acid_sites = [s for s in result.sites if s.kind == "acid"]
    assert acid_sites
    closest = min(acid_sites, key=lambda s: abs(s.pka - expected_pka))
    assert closest.pka == pytest.approx(expected_pka, abs=1.5)


@pytest.mark.parametrize(
    "smiles,expected_pkah",
    [
        (ANILINE, 4.6),
        (PYRIDINE, 5.2),
        (IMIDAZOLE, 7.0),
    ],
)
def test_base_pka_matches_literature(predictor, smiles, expected_pkah):
    result = predictor.predict(smiles, ph=7.4, max_protomers=5)
    assert result is not None
    base_sites = [s for s in result.sites if s.kind == "base"]
    assert base_sites
    closest = min(base_sites, key=lambda s: abs(s.pka - expected_pkah))
    assert closest.pka == pytest.approx(expected_pkah, abs=1.5)


def test_caffeine_has_no_acid_sites(predictor):
    """Negative control: caffeine has no acidic protons. The fix must not
    spuriously invent acid sites where none exist."""
    result = predictor.predict(CAFFEINE, ph=7.4, max_protomers=5)
    assert result is not None
    assert [s for s in result.sites if s.kind == "acid"] == []
    assert result.microspecies[0].charge == 0
    assert result.microspecies[0].population > 0.99


def test_symmetric_amine_sites_are_deduped(predictor):
    """Regression test for the duplicate-acid-site bug: aniline's -NH2 has
    two equivalent H's: before the dedup fix this produced two identical
    acid sites at the same atom_idx."""
    result = predictor.predict(ANILINE, ph=7.4, max_protomers=8)
    assert result is not None
    keys = [(s.kind, s.atom_idx) for s in result.sites]
    assert len(keys) == len(set(keys)), f"duplicate (kind, atom_idx) sites: {keys}"


def test_sulfanilamide_sites_are_deduped(predictor):
    """Sulfanilamide has two -NH2-type groups (aniline + sulfonamide), so
    it's a stronger stress test for the dedup than aniline alone."""
    result = predictor.predict(SULFANILAMIDE, ph=7.4, max_protomers=8)
    assert result is not None
    keys = [(s.kind, s.atom_idx) for s in result.sites]
    assert len(keys) == len(set(keys)), f"duplicate (kind, atom_idx) sites: {keys}"
    assert len(result.sites) == 4  # 2 acid (aniline N, sulfonamide N) + 2 base


def test_microspecies_populations_sum_to_one(predictor):
    for smiles in [ACETAMINOPHEN, SULFANILAMIDE, ANILINE, IMIDAZOLE]:
        result = predictor.predict(smiles, ph=7.4, max_protomers=8)
        assert result is not None
        total = sum(m.population for m in result.microspecies)
        assert total == pytest.approx(1.0, abs=1e-6)


def test_imidazole_protonation_fraction_matches_henderson_hasselbalch(predictor):
    """Imidazole's own pKaH (~7.0) sits close to physiological pH, so this
    also checks the independent-site population math, not just site
    detection: fraction protonated should follow 1 / (1 + 10^(pH - pKa))."""
    result = predictor.predict(IMIDAZOLE, ph=7.4, max_protomers=5)
    assert result is not None
    base_site = min((s for s in result.sites if s.kind == "base"), key=lambda s: s.pka)
    expected_frac = 1.0 / (1.0 + 10.0 ** (7.4 - base_site.pka))
    cation_pop = sum(m.population for m in result.microspecies if m.charge == 1)
    assert cation_pop == pytest.approx(expected_frac, abs=0.05)
