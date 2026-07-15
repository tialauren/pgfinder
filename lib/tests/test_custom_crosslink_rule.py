"""Test the user-supplied custom regex crosslink rule (for species not in the built-in list)."""

from decimal import Decimal

import pandas as pd
import pytest

from pgfinder import COLUMNS
from pgfinder.accurate_dimer_builder import AccurateDimerBuilder, build_dimers_for_species
from pgfinder.accurate_species_rules import CustomCrosslinkRule, match_custom_acceptor, match_custom_donor

DEFAULT_COLUMNS = COLUMNS["pgfinder"]

ECOLI_LIKE_RULE = CustomCrosslinkRule(
    label="Custom 4-3",
    donor_pattern=r"^..JAA$",
    acceptor_pattern=r"^..J",
    loses_terminal_ala=True,
)


def test_match_custom_donor_and_acceptor():
    assert match_custom_donor("gm-AEJAA|1", ECOLI_LIKE_RULE) is True
    assert match_custom_acceptor("gm-AEJAA|1", ECOLI_LIKE_RULE) is True
    # A tripeptide doesn't end in JAA, so it isn't a valid donor...
    assert match_custom_donor("gm-AEJ|1", ECOLI_LIKE_RULE) is False
    # ...but it still has J at position 3, so it's a valid acceptor.
    assert match_custom_acceptor("gm-AEJ|1", ECOLI_LIKE_RULE) is True
    # A Lys-type stem (K, not J) matches neither.
    assert match_custom_donor("gm-AQKAA|1", ECOLI_LIKE_RULE) is False
    assert match_custom_acceptor("gm-AQKAA|1", ECOLI_LIKE_RULE) is False


def test_match_custom_rule_handles_unparseable_structure():
    assert match_custom_donor("not a valid structure", ECOLI_LIKE_RULE) is False
    assert match_custom_acceptor("not a valid structure", ECOLI_LIKE_RULE) is False


def test_match_custom_rule_handles_invalid_regex():
    bad_rule = CustomCrosslinkRule(
        label="Bad", donor_pattern="[unterminated", acceptor_pattern="[unterminated", loses_terminal_ala=False
    )
    assert match_custom_donor("gm-AEJAA|1", bad_rule) is False
    assert match_custom_acceptor("gm-AEJAA|1", bad_rule) is False


def test_dimer_builder_requires_exactly_one_of_species_code_or_custom_rule():
    with pytest.raises(ValueError):
        AccurateDimerBuilder()
    with pytest.raises(ValueError):
        AccurateDimerBuilder(species_code="ecoli", custom_rule=ECOLI_LIKE_RULE)


def test_custom_rule_dimer_matches_equivalent_built_in_species_rule():
    """A custom rule mimicking E. coli's 4-3 mechanism should produce the same
    trimming and mass as the built-in 'ecoli' species rule."""
    donor_mass = Decimal("1012.444817")
    acceptor_mass = Decimal("1012.444817")
    acceptors_df = pd.DataFrame([{"structure": "gm-AEJAA|1", "mass": str(acceptor_mass)}])

    builder = AccurateDimerBuilder(custom_rule=ECOLI_LIKE_RULE)
    dimers = builder._generate_dimers(
        donor_structure="gm-AEJAA|1",
        donor_mass=donor_mass,
        acceptors_df=acceptors_df,
        crosslink_type="Custom 4-3",
        columns=DEFAULT_COLUMNS,
    )

    assert len(dimers) == 1
    dimer = dimers[0]
    assert dimer["Inferred structure"] == "gm-AEJA=gm-AEJAA[Custom 4-3]|2"

    water_mass = Decimal("18.0106")
    d_ala_mass = Decimal("89.047679")
    expected_mass = float(donor_mass + acceptor_mass - water_mass - d_ala_mass)
    assert dimer["Theo (Da)"] == pytest.approx(expected_mass, abs=1e-6)


def test_custom_rule_without_lost_residue_keeps_full_donor():
    rule = CustomCrosslinkRule(
        label="Custom 3-3", donor_pattern=r"^..J.$", acceptor_pattern=r"^..J", loses_terminal_ala=False
    )
    donor_mass = Decimal("941.407703")
    acceptors_df = pd.DataFrame([{"structure": "gm-AEJA|1", "mass": str(donor_mass)}])

    builder = AccurateDimerBuilder(custom_rule=rule)
    dimers = builder._generate_dimers(
        donor_structure="gm-AEJA|1",
        donor_mass=donor_mass,
        acceptors_df=acceptors_df,
        crosslink_type="Custom 3-3",
        columns=DEFAULT_COLUMNS,
    )

    # Custom labels aren't recognized by the '~'-for-3-3 joiner convention
    # (that's keyed on the literal string '3-3'), so they fall back to '='.
    assert dimers[0]["Inferred structure"] == "gm-AEJA=gm-AEJA[Custom 3-3]|2"
    expected_mass = float(donor_mass * 2 - Decimal("18.0106"))
    assert dimers[0]["Theo (Da)"] == pytest.approx(expected_mass, abs=1e-6)


@pytest.mark.parametrize(
    "structure,expected",
    [
        ("gm-AQK[GGGGG]AA|1", True),  # exactly 5 glycines - matches
        ("gm-AQK[GGG]AA|1", False),  # 3 glycines - outside [5,5]
        ("gm-AQKAA|1", False),  # no bridge at all
    ],
)
def test_glycine_bridge_acceptor_matching(structure, expected):
    rule = CustomCrosslinkRule(
        label="Custom bridge",
        donor_pattern=r"^..KAA$",
        acceptor_pattern=r"^..K",
        loses_terminal_ala=True,
        acceptor_bridge_type="glycine",
        acceptor_min_glycine_bridge=5,
        acceptor_max_glycine_bridge=5,
    )
    assert match_custom_acceptor(structure, rule) is expected


def test_glycine_bridge_acceptor_matching_permissive_range():
    rule = CustomCrosslinkRule(
        label="Custom bridge",
        donor_pattern=r"^..KAA$",
        acceptor_pattern=r"^..K",
        loses_terminal_ala=True,
        acceptor_bridge_type="glycine",
        acceptor_min_glycine_bridge=1,
        acceptor_max_glycine_bridge=5,
    )
    assert match_custom_acceptor("gm-AQK[GGG]AA|1", rule) is True
    assert match_custom_acceptor("gm-AQK[GGGGGG]AA|1", rule) is False  # 6 glycines - outside [1,5]


@pytest.mark.parametrize(
    "structure,expected",
    [
        ("gm-AQK[D]AA|1", True),
        ("gm-AQKAA|1", False),
        ("gm-AQK[GGGGG]AA|1", False),  # glycine bridge, not D-Asp
    ],
)
def test_dasp_bridge_acceptor_matching(structure, expected):
    rule = CustomCrosslinkRule(
        label="Custom dasp",
        donor_pattern=r"^..KAA$",
        acceptor_pattern=r"^..K",
        loses_terminal_ala=True,
        acceptor_bridge_type="dasp",
    )
    assert match_custom_acceptor(structure, rule) is expected


def test_no_bridge_type_ignores_bridge_info():
    """acceptor_bridge_type='none' (the default) should match regardless of
    whether the structure happens to have a bridge - same behavior as before
    bridge-awareness was added."""
    rule = CustomCrosslinkRule(
        label="Custom no-bridge-check", donor_pattern=r"^..KAA$", acceptor_pattern=r"^..K", loses_terminal_ala=True
    )
    assert match_custom_acceptor("gm-AQKAA|1", rule) is True
    assert match_custom_acceptor("gm-AQK[GGGGG]AA|1", rule) is True


def test_build_dimers_for_species_with_custom_rule_end_to_end():
    columns = DEFAULT_COLUMNS
    matched_monomers = pd.DataFrame(
        [
            {
                columns["inferred"]["structure"]: "gm-AEJAA|1",
                columns["inferred"]["mass"]: 1012.444817,
                columns["input"]["intensity"]: 5000.0,
            }
        ]
    )
    theo_db = pd.DataFrame([{columns["inferred"]["structure"]: "gm-AEJAA|1", columns["inferred"]["mass"]: 1012.444817}])

    result, donors_used = build_dimers_for_species(matched_monomers, theo_db, custom_rule=ECOLI_LIKE_RULE)

    assert not result.empty
    assert (result["Inferred structure"] == "gm-AEJA=gm-AEJAA[Custom 4-3]|2").any()
    assert (donors_used["Crosslink Type"] == "Custom 4-3").any()


def test_custom_bridge_rule_matches_built_in_saureus_rule():
    """A custom rule mimicking S. aureus's pentaglycine-bridge mechanism should
    accept/reject the same structures as the built-in 'saureus' species rule."""
    from pgfinder.accurate_species_rules import AccurateSpeciesRules

    rule = CustomCrosslinkRule(
        label="Custom saureus",
        donor_pattern=r"^..KAA$",
        acceptor_pattern=r"^..K",
        loses_terminal_ala=True,
        acceptor_bridge_type="glycine",
        acceptor_min_glycine_bridge=5,
        acceptor_max_glycine_bridge=5,
    )

    structures = [
        "gm-AQKAA|1",
        "gm-AQK[GGGGG]AA|1",
        "gm-AQK[GGG]AA|1",
        "gm-AEJAA|1",
    ]
    for structure in structures:
        built_in = AccurateSpeciesRules.is_valid_acceptor(structure, "saureus", "4-3-bridge", strict_mode=True)
        custom = match_custom_acceptor(structure, rule)
        assert built_in == custom, f"mismatch for {structure}: built_in={built_in}, custom={custom}"
