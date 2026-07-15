"""Test the accurate dimer builder.

Covers two things that are easy to silently break:
- mass accuracy (the lost D-Ala must use the monoisotopic mass convention
  that every other "Theo (Da)" value in pgfinder uses)
- that column names are read from the `columns` mapping rather than
  hardcoded literals, so callers using a non-default column config still work
"""

from decimal import Decimal

import pandas as pd
import pytest

from pgfinder import COLUMNS
from pgfinder.accurate_dimer_builder import AccurateDimerBuilder, build_dimers_for_species

DEFAULT_COLUMNS = COLUMNS["pgfinder"]

WATER_MASS = Decimal("18.0106")
D_ALA_MONOISOTOPIC_MASS = Decimal("89.047679")


def test_generate_dimers_trims_donor_and_uses_monoisotopic_mass():
    """4-3 crosslinks must trim the donor to a tetrapeptide and lose one
    monoisotopic D-Ala, not the average-mass value."""
    builder = AccurateDimerBuilder("ecoli")
    donor_mass = Decimal("1012.444817")
    acceptor_mass = Decimal("1012.444817")
    acceptors_df = pd.DataFrame([{"structure": "gm-AEJAA|1", "mass": str(acceptor_mass)}])

    dimers = builder._generate_dimers(
        donor_structure="gm-AEJAA|1",
        donor_mass=donor_mass,
        acceptors_df=acceptors_df,
        crosslink_type="4-3",
        columns=DEFAULT_COLUMNS,
    )

    assert len(dimers) == 1
    dimer = dimers[0]
    # Donor is trimmed from the pentapeptide (AEJAA) to the tetrapeptide (AEJA).
    assert dimer["Inferred structure"] == "gm-AEJA=gm-AEJAA[4-3]|2"

    expected_mass = float(donor_mass + acceptor_mass - WATER_MASS - D_ALA_MONOISOTOPIC_MASS)
    assert dimer["Theo (Da)"] == pytest.approx(expected_mass, abs=1e-6)


def test_generate_dimers_3_3_crosslink_trims_donor_to_tripeptide():
    """3-3 crosslinks: an L,D-transpeptidase cleaves the tetrapeptide donor's
    mDAP3-D-Ala4 bond, releasing free D-Ala4 - same lost-residue mechanism as
    4-3, just one position earlier. The acceptor is unconstrained in length."""
    builder = AccurateDimerBuilder("ecoli")
    donor_mass = Decimal("941.407703")
    acceptor_mass = Decimal("941.407703")
    acceptors_df = pd.DataFrame([{"structure": "gm-AEJA|1", "mass": str(acceptor_mass)}])

    dimers = builder._generate_dimers(
        donor_structure="gm-AEJA|1",
        donor_mass=donor_mass,
        acceptors_df=acceptors_df,
        crosslink_type="3-3",
        columns=DEFAULT_COLUMNS,
    )

    assert len(dimers) == 1
    dimer = dimers[0]
    # Donor is trimmed from the tetrapeptide (AEJA) to the tripeptide (AEJ);
    # the acceptor keeps its full length (AEJA).
    assert dimer["Inferred structure"] == "gm-AEJ=gm-AEJA[3-3]|2"

    expected_mass = float(donor_mass + acceptor_mass - WATER_MASS - D_ALA_MONOISOTOPIC_MASS)
    assert dimer["Theo (Da)"] == pytest.approx(expected_mass, abs=1e-6)


def test_build_dimers_respects_custom_columns_mapping():
    """Column names must come from the `columns` argument, not be hardcoded."""
    custom_columns = {
        "input": {"intensity": "MyIntensity"},
        "inferred": {"structure": "MyStructure", "mass": "MyMass"},
    }
    matched_monomers = pd.DataFrame([{"MyStructure": "gm-AEJAA|1", "MyMass": 1012.444817, "MyIntensity": 5000.0}])
    theo_db = pd.DataFrame([{"MyStructure": "gm-AEJAA|1", "MyMass": 1012.444817}])

    builder = AccurateDimerBuilder("ecoli")
    result, donors_used = builder.build_dimers_from_detected_monomers(matched_monomers, theo_db, columns=custom_columns)

    assert list(result.columns) == ["MyStructure", "MyMass"]
    assert len(result) >= 1
    assert "MyStructure" in donors_used.columns


def test_build_dimers_for_species_ecoli_real_structures():
    """End-to-end smoke test using real E. coli structure notation (J = mDAP)."""
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

    result, donors_used = build_dimers_for_species(matched_monomers, theo_db, "ecoli")

    assert not result.empty
    assert (result["Inferred structure"] == "gm-AEJA=gm-AEJAA[4-3]|2").any()
    assert len(donors_used) == 1
    assert donors_used.iloc[0]["Inferred structure"] == "gm-AEJAA|1"
    assert donors_used.iloc[0]["Percent of eligible total"] == pytest.approx(100.0)


def test_generate_dimers_trims_modified_donor_and_preserves_tag():
    """A donor carrying a modification tag (e.g. Anhydro-MurNAc) must still be
    trimmed correctly, with the tag re-attached after the trimmed stem rather
    than lost or left attached to the untrimmed residue."""
    builder = AccurateDimerBuilder("ecoli")
    donor_mass = Decimal("992.418617")  # gm-AEJAA mass minus Anhydro-MurNAc's -20.0262
    acceptor_mass = Decimal("1012.444817")
    acceptors_df = pd.DataFrame([{"structure": "gm-AEJAA|1", "mass": str(acceptor_mass)}])

    dimers = builder._generate_dimers(
        donor_structure="gm-AEJAA (Anh)|1",
        donor_mass=donor_mass,
        acceptors_df=acceptors_df,
        crosslink_type="4-3",
        columns=DEFAULT_COLUMNS,
    )

    assert len(dimers) == 1
    assert dimers[0]["Inferred structure"] == "gm-AEJA (Anh)=gm-AEJAA[4-3]|2"

    expected_mass = float(donor_mass + acceptor_mass - WATER_MASS - D_ALA_MONOISOTOPIC_MASS)
    assert dimers[0]["Theo (Da)"] == pytest.approx(expected_mass, abs=1e-6)


def test_build_dimers_for_species_prefers_more_abundant_modified_donor():
    """Before `split_modification_tag` existed, a modified monomer failed to
    parse entirely and was silently invisible to donor selection - meaning a
    less-abundant unmodified monomer would be picked instead of the true most
    abundant donor. This is the regression test for that fix."""
    columns = DEFAULT_COLUMNS
    matched_monomers = pd.DataFrame(
        [
            {
                columns["inferred"]["structure"]: "gm-AEJAA (Anh)|1",
                columns["inferred"]["mass"]: 992.418617,
                columns["input"]["intensity"]: 9000.0,
            },
            {
                columns["inferred"]["structure"]: "gm-AEJAA|1",
                columns["inferred"]["mass"]: 1012.444817,
                columns["input"]["intensity"]: 1000.0,
            },
        ]
    )
    theo_db = pd.DataFrame([{columns["inferred"]["structure"]: "gm-AEJAA|1", columns["inferred"]["mass"]: 1012.444817}])

    result, donors_used = build_dimers_for_species(matched_monomers, theo_db, "ecoli")

    assert not result.empty
    assert (result["Inferred structure"] == "gm-AEJA (Anh)=gm-AEJAA[4-3]|2").any()
    assert (donors_used["Inferred structure"] == "gm-AEJAA (Anh)|1").any()


# ---------------------------------------------------------------------------
# Donor abundance threshold tests
# ---------------------------------------------------------------------------


def _multi_donor_monomers(columns=DEFAULT_COLUMNS):
    """Three monomers for threshold testing (E. coli 4-3 context):
    - gm-AEJAA|1         70 000 intensity  (70 % of the 4-3 eligible pool)
    - gm-AEJAA (Anh)|1   30 000 intensity  (30 % of the 4-3 eligible pool)
    - gm-AEJA|1         900 000 intensity  (ineligible for 4-3 — tetrapeptide;
                                             eligible for 3-3 as the sole donor)
    """
    return pd.DataFrame(
        [
            {
                columns["inferred"]["structure"]: "gm-AEJAA|1",
                columns["inferred"]["mass"]: 1012.444817,
                columns["input"]["intensity"]: 70_000.0,
            },
            {
                columns["inferred"]["structure"]: "gm-AEJAA (Anh)|1",
                columns["inferred"]["mass"]: 992.418617,
                columns["input"]["intensity"]: 30_000.0,
            },
            {
                columns["inferred"]["structure"]: "gm-AEJA|1",
                columns["inferred"]["mass"]: 941.407703,
                columns["input"]["intensity"]: 900_000.0,
            },
        ]
    )


def test_threshold_low_always_includes_at_least_one_donor():
    """Even at threshold 0.0, the most abundant eligible donor is always
    included — this recovers the single-donor behaviour as the limit case."""
    builder = AccurateDimerBuilder("ecoli", donor_abundance_threshold=0.0)
    donors = builder._find_eligible_donors_within_threshold(_multi_donor_monomers(), "4-3", DEFAULT_COLUMNS)

    assert len(donors) == 1
    assert donors[0]["structure"] == "gm-AEJAA|1"


def test_threshold_100_includes_all_eligible_donors():
    """At threshold=1.0 every structurally eligible donor is returned, even
    low-abundance ones that would normally not be worth including."""
    builder = AccurateDimerBuilder("ecoli", donor_abundance_threshold=1.0)
    donors = builder._find_eligible_donors_within_threshold(_multi_donor_monomers(), "4-3", DEFAULT_COLUMNS)

    assert len(donors) == 2
    assert {d["structure"] for d in donors} == {"gm-AEJAA|1", "gm-AEJAA (Anh)|1"}


def test_ineligible_monomers_excluded_from_eligible_pool_denominator():
    """A very abundant but structurally ineligible monomer (gm-AEJA|1 is a
    tetrapeptide and cannot be a 4-3 donor) must not appear in results and must
    not count toward the eligible pool's total intensity. If it did, percentages
    would be distorted and the threshold logic would misbehave."""
    builder = AccurateDimerBuilder("ecoli", donor_abundance_threshold=1.0)
    donors = builder._find_eligible_donors_within_threshold(_multi_donor_monomers(), "4-3", DEFAULT_COLUMNS)

    structures = {d["structure"] for d in donors}
    assert "gm-AEJA|1" not in structures

    # Eligible pool = 70 000 + 30 000 = 100 000 (not 1 000 000)
    top = next(d for d in donors if d["structure"] == "gm-AEJAA|1")
    second = next(d for d in donors if "Anh" in d["structure"])
    assert top["percent_of_eligible_total"] == pytest.approx(70.0)
    assert second["percent_of_eligible_total"] == pytest.approx(30.0)
    assert top["cumulative_percent"] == pytest.approx(70.0)
    assert second["cumulative_percent"] == pytest.approx(100.0)


def test_donors_used_df_percent_math_is_correct():
    """Percent of eligible total and cumulative percent in donors_used_df must
    be computed relative to the eligible pool only and accumulate correctly."""
    theo_db = pd.DataFrame(
        [{DEFAULT_COLUMNS["inferred"]["structure"]: "gm-AEJAA|1", DEFAULT_COLUMNS["inferred"]["mass"]: 1012.444817}]
    )
    _, donors_used = build_dimers_for_species(_multi_donor_monomers(), theo_db, "ecoli", donor_abundance_threshold=1.0)

    four_three = donors_used[donors_used["Crosslink Type"] == "4-3 crosslink"].copy()
    assert len(four_three) == 2

    four_three.sort_values(by=DEFAULT_COLUMNS["input"]["intensity"], ascending=False, inplace=True)
    top = four_three.iloc[0]
    second = four_three.iloc[1]

    assert top["Inferred structure"] == "gm-AEJAA|1"
    assert top["Percent of eligible total"] == pytest.approx(70.0)
    assert top["Cumulative percent"] == pytest.approx(70.0)
    assert second["Percent of eligible total"] == pytest.approx(30.0)
    assert second["Cumulative percent"] == pytest.approx(100.0)


def test_multiple_donors_each_generate_their_own_dimers():
    """When the threshold covers multiple donors, dimers from every selected
    donor must appear in the output — not just from the most abundant one."""
    theo_db = pd.DataFrame(
        [{DEFAULT_COLUMNS["inferred"]["structure"]: "gm-AEJAA|1", DEFAULT_COLUMNS["inferred"]["mass"]: 1012.444817}]
    )
    result, donors_used = build_dimers_for_species(
        _multi_donor_monomers(), theo_db, "ecoli", donor_abundance_threshold=1.0
    )

    four_three_donors = donors_used[donors_used["Crosslink Type"] == "4-3 crosslink"]
    assert len(four_three_donors) == 2

    structures = set(result["Inferred structure"])
    # 4-3 dimer from unmodified donor (gm-AEJAA trimmed to gm-AEJA)
    assert "gm-AEJA=gm-AEJAA[4-3]|2" in structures
    # 4-3 dimer from modified donor (gm-AEJAA (Anh) trimmed to gm-AEJA (Anh))
    assert "gm-AEJA (Anh)=gm-AEJAA[4-3]|2" in structures
    # 3-3 dimer from gm-AEJA|1 (only eligible 3-3 donor, trimmed to gm-AEJ)
    assert "gm-AEJ=gm-AEJAA[3-3]|2" in structures
