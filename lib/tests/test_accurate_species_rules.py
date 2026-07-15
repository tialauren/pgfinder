"""Test biochemically accurate species/crosslink rules.

These are regression tests for structure strings that actually appear in
pgfinder's own theoretical mass databases (lib/pgfinder/masses/*.csv), not
contrived examples - they exercise the real "J" (mDAP) / "K" (Lysine) and
"g(-Ac)m-" glycan-prefix notation used by the project.
"""

import pytest

from pgfinder.accurate_species_rules import AccurateSpeciesRules, parse_peptide_structure, split_modification_tag


def test_parse_peptide_structure_ecoli_pentapeptide():
    """Real E. coli pentapeptide monomer: mDAP is denoted 'J' at position 3."""
    info = parse_peptide_structure("gm-AEJAA|1")
    assert info.sequence == ["A", "E", "J", "A", "A"]
    assert info.length == 5
    assert info.has_mdap_pos3 is True
    assert info.has_lys_pos3 is False


def test_parse_peptide_structure_lys_type_pentapeptide():
    """Lys-type pentapeptide (e.g. S. aureus): Lysine is denoted 'K' at position 3."""
    info = parse_peptide_structure("gm-AQKAA|1")
    assert info.sequence == ["A", "Q", "K", "A", "A"]
    assert info.has_mdap_pos3 is False
    assert info.has_lys_pos3 is True


def test_parse_peptide_structure_deacetylated_prefix():
    """C. diff structures use the 'g(-Ac)m-' de-N-acetylated glycan prefix."""
    info = parse_peptide_structure("g(-Ac)m-AEJAA|1")
    assert info is not None
    assert info.sequence == ["A", "E", "J", "A", "A"]
    assert info.has_mdap_pos3 is True


def test_parse_peptide_structure_with_glycine_bridge():
    info = parse_peptide_structure("gm-AQK[GGGGG]AA|1")
    assert info.sequence == ["A", "Q", "K", "A", "A"]
    assert info.glycine_bridge_length == 5


@pytest.mark.parametrize(
    "structure,species,crosslink_type,expected",
    [
        ("gm-AEJAA|1", "ecoli", "4-3", True),
        ("g(-Ac)m-AEJAA|1", "cdiff", "4-3", True),
        ("gm-AQKAA|1", "saureus", "4-3-bridge", True),
        # A Lys-type stem is not a valid E. coli (DAP-type) donor.
        ("gm-AQKAA|1", "ecoli", "4-3", False),
        # A mDAP-type stem is not a valid S. aureus (Lys-type) donor.
        ("gm-AEJAA|1", "saureus", "4-3-bridge", False),
    ],
)
def test_is_valid_donor(structure, species, crosslink_type, expected):
    assert AccurateSpeciesRules.is_valid_donor(structure, species, crosslink_type) is expected


def test_is_valid_acceptor_saureus_strict_requires_five_glycines():
    acceptor_full_bridge = "gm-AQK[GGGGG]AA|1"
    acceptor_short_bridge = "gm-AQK[GGG]AA|1"

    assert AccurateSpeciesRules.is_valid_acceptor(acceptor_full_bridge, "saureus", "4-3-bridge", strict_mode=True)
    assert not AccurateSpeciesRules.is_valid_acceptor(acceptor_short_bridge, "saureus", "4-3-bridge", strict_mode=True)
    # Permissive mode allows shorter bridges for novel-crosslink discovery.
    assert AccurateSpeciesRules.is_valid_acceptor(acceptor_short_bridge, "saureus", "4-3-bridge", strict_mode=False)


@pytest.mark.parametrize(
    "structure_base,expected_stem,expected_tag",
    [
        ("gm-AEJAA (Anh)", "gm-AEJAA", " (Anh)"),
        ("gm-AEJAA (-Ac, Anh)", "gm-AEJAA", " (-Ac, Anh)"),
        ("gm-AEJAA", "gm-AEJAA", ""),
        # The C. diff de-acetylated glycan prefix isn't a trailing tag and must survive untouched.
        ("g(-Ac)m-AEJAA", "g(-Ac)m-AEJAA", ""),
    ],
)
def test_split_modification_tag(structure_base, expected_stem, expected_tag):
    assert split_modification_tag(structure_base) == (expected_stem, expected_tag)


@pytest.mark.parametrize(
    "structure",
    [
        "gm-AEJAA (Anh)|1",
        "gm-AEJAA (Am)|1",
        "gm-AEJAA (-Ac)|1",
        "gm-AEJAA (-Ac, Anh)|1",
        "g(-Ac)m-AEJAA (Anh)|1",
    ],
)
def test_parse_peptide_structure_strips_modification_tags(structure):
    """A modification tag (appended by `matching.modification_generator`) must
    not prevent the stem from being parsed - modifications act on the glycan
    or side chains, not on stem length/position-3 identity."""
    info = parse_peptide_structure(structure)
    assert info is not None
    assert info.sequence == ["A", "E", "J", "A", "A"]
    assert info.has_mdap_pos3 is True


def test_is_valid_donor_recognizes_modified_structure():
    """A modified monomer must still be classified as a valid donor - this was
    silently broken before `split_modification_tag` existed (the tag broke
    parsing entirely, so any modified monomer was invisible to donor/acceptor
    selection, not incorrectly priced)."""
    assert AccurateSpeciesRules.is_valid_donor("gm-AEJAA (Anh)|1", "ecoli", "4-3") is True
