"""Test cumulative-abundance donor selection.

`AccurateDimerBuilder` used to pick exactly one donor per crosslink type (the
single most abundant structurally-eligible monomer). It now selects however
many of the most abundant *eligible* donors are needed to cover
`donor_abundance_threshold` of the eligible pool's cumulative intensity -
these tests cover that algorithm directly, independent of the dimers it
goes on to build.
"""

import pandas as pd
import pytest

from pgfinder import COLUMNS
from pgfinder.accurate_dimer_builder import AccurateDimerBuilder

DEFAULT_COLUMNS = COLUMNS["pgfinder"]

# Four structures that are all valid E. coli 4-3 donors (length 5, mDAP at
# position 3, ending "AA"), at decreasing intensity: 50%, 30%, 15%, 5% of the
# combined total. Mixed in is one high-intensity Lys-type monomer that is
# *not* eligible for ecoli's mDAP-type donor rule.
FOUR_THREE_DONORS = pd.DataFrame(
    [
        {"Inferred structure": "gm-AEJAA|1", "Theo (Da)": 1012.0, "Intensity": 5000.0},
        {"Inferred structure": "gm-AQJAA|1", "Theo (Da)": 1020.0, "Intensity": 3000.0},
        {"Inferred structure": "gm-AFJAA|1", "Theo (Da)": 1030.0, "Intensity": 1500.0},
        {"Inferred structure": "gm-AGJAA|1", "Theo (Da)": 1040.0, "Intensity": 500.0},
        # Ineligible (Lys-type) but more abundant than every eligible donor -
        # must not shrink the eligible pool's denominator.
        {"Inferred structure": "gm-AQKAA|1", "Theo (Da)": 1050.0, "Intensity": 50000.0},
    ]
)


def test_low_threshold_selects_single_most_abundant_donor():
    """A low threshold recovers the old single-donor behavior as a special
    case - at least the most abundant eligible donor is always included."""
    builder = AccurateDimerBuilder("ecoli", donor_abundance_threshold=0.01)
    donors = builder._find_eligible_donors_within_threshold(FOUR_THREE_DONORS, "4-3", DEFAULT_COLUMNS)

    assert len(donors) == 1
    assert donors[0]["structure"] == "gm-AEJAA|1"
    assert donors[0]["percent_of_eligible_total"] == pytest.approx(50.0)


def test_threshold_selects_multiple_donors_until_cumulative_covers_it():
    builder = AccurateDimerBuilder("ecoli", donor_abundance_threshold=0.9)
    donors = builder._find_eligible_donors_within_threshold(FOUR_THREE_DONORS, "4-3", DEFAULT_COLUMNS)

    # 50% + 30% + 15% = 95%, crossing the 90% threshold on the third donor.
    assert [d["structure"] for d in donors] == ["gm-AEJAA|1", "gm-AQJAA|1", "gm-AFJAA|1"]
    assert donors[-1]["cumulative_percent"] == pytest.approx(95.0)


def test_full_threshold_selects_every_eligible_donor():
    builder = AccurateDimerBuilder("ecoli", donor_abundance_threshold=1.0)
    donors = builder._find_eligible_donors_within_threshold(FOUR_THREE_DONORS, "4-3", DEFAULT_COLUMNS)

    assert len(donors) == 4
    assert donors[-1]["cumulative_percent"] == pytest.approx(100.0)


def test_ineligible_abundant_monomer_excluded_from_pool_and_denominator():
    """The Lys-type monomer (intensity 50000, by far the largest) must not be
    selected as a 4-3 donor, and must not dilute the percentages of the
    actually-eligible donors - percentages are computed over the eligible
    pool only (5000+3000+1500+500 = 10000), not the full 60000."""
    builder = AccurateDimerBuilder("ecoli", donor_abundance_threshold=1.0)
    donors = builder._find_eligible_donors_within_threshold(FOUR_THREE_DONORS, "4-3", DEFAULT_COLUMNS)

    assert all(d["structure"] != "gm-AQKAA|1" for d in donors)
    assert donors[0]["percent_of_eligible_total"] == pytest.approx(50.0)  # not 5000/60000


def test_no_eligible_donors_returns_empty_list():
    mdap_only = FOUR_THREE_DONORS[FOUR_THREE_DONORS["Inferred structure"] != "gm-AQKAA|1"]
    builder = AccurateDimerBuilder("saureus", donor_abundance_threshold=0.9)
    # None of these are Lys-type-with-pentaglycine-bridge donors (all mDAP-type).
    donors = builder._find_eligible_donors_within_threshold(mdap_only, "4-3-bridge", DEFAULT_COLUMNS)

    assert donors == []


def test_duplicate_structure_rows_are_aggregated_not_treated_as_separate_donors():
    """The same structure can appear as multiple rows (e.g. matched at
    different retention times/charge states) - these must be combined into a
    single donor with summed intensity, not treated as two distinct
    candidates (which would otherwise generate duplicate dimers downstream)."""
    duplicated = pd.DataFrame(
        [
            {"Inferred structure": "gm-AEJAA|1", "Theo (Da)": 1012.0, "Intensity": 3000.0},
            {"Inferred structure": "gm-AEJAA|1", "Theo (Da)": 1012.0, "Intensity": 2000.0},
            {"Inferred structure": "gm-AQJAA|1", "Theo (Da)": 1020.0, "Intensity": 1000.0},
        ]
    )
    builder = AccurateDimerBuilder("ecoli", donor_abundance_threshold=1.0)
    donors = builder._find_eligible_donors_within_threshold(duplicated, "4-3", DEFAULT_COLUMNS)

    assert len(donors) == 2
    assert donors[0]["structure"] == "gm-AEJAA|1"
    assert donors[0]["intensity"] == pytest.approx(5000.0)  # 3000 + 2000 combined
    assert donors[0]["percent_of_eligible_total"] == pytest.approx(5000 / 6000 * 100)
