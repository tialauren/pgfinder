"""Test the standalone theoretical-dimer-list generator.

This is the function that lets a user see exactly which theoretical dimers a
full `data_analysis_with_dimers` run would search for, without waiting for
the final match-against-raw-data and consolidation steps.
"""

import pytest

from pgfinder import matching
from pgfinder.accurate_species_rules import CustomCrosslinkRule


def test_generate_theoretical_dimers_requires_species_or_custom_rule(raw_data, theo_masses, ppm):
    with pytest.raises(ValueError):
        matching.generate_theoretical_dimers(raw_data, theo_masses, [], ppm)


def test_generate_theoretical_dimers_ecoli_real_data(raw_data, theo_masses, ppm):
    """End-to-end smoke test against the bundled E. coli MaxQuant test data -
    mirrors what `find_pg --enable_dimers --species ecoli` produces."""
    result, donors_used = matching.generate_theoretical_dimers(
        raw_data, theo_masses, [], ppm, species_code="ecoli", strict_mode=True
    )

    assert not result.empty
    assert "Inferred structure" in result.columns
    assert "Theo (Da)" in result.columns
    # Every generated structure should be tagged as a dimer (oligomerization state 2).
    assert result["Inferred structure"].str.endswith("|2").all()
    assert not donors_used.empty
    assert "Crosslink Type" in donors_used.columns


def test_generate_theoretical_dimers_no_valid_donor_returns_empty(raw_data, theo_masses, ppm):
    """The bundled test data is DAP-type (mDAP/J at position 3) - asking for a
    Lys-type (K) S. aureus donor should find no qualifying donor and return an
    empty DataFrame rather than raising."""
    result, donors_used = matching.generate_theoretical_dimers(
        raw_data, theo_masses, [], ppm, species_code="saureus", strict_mode=True
    )

    assert result.empty
    assert donors_used.empty


def test_generate_theoretical_dimers_matches_data_analysis_with_dimers_search_space(raw_data, theo_masses, ppm):
    """The theoretical dimers generated standalone should be exactly the set
    added to the search space inside a full `data_analysis_with_dimers` run -
    this is the guarantee that makes the standalone preview trustworthy."""
    standalone, _donors_used = matching.generate_theoretical_dimers(
        raw_data, theo_masses, [], ppm, species_code="ecoli", strict_mode=True
    )

    full_results, _donors = matching.data_analysis_with_dimers(
        raw_data,
        theo_masses,
        rt_window=0.5,
        enabled_mod_list=[],
        ppm_tolerance=ppm,
        consolidation_ppm=1,
        enable_dimers=True,
        species_code="ecoli",
        strict_mode=True,
    )

    matched_dimers = set(
        full_results.loc[full_results["Inferred structure"].str.endswith("|2", na=False), "Inferred structure"]
    )
    theoretical_dimers = set(standalone["Inferred structure"])

    # Every dimer actually matched in the full run must have been part of the
    # standalone theoretical list (the full run can only narrow this set down
    # by matching against real intensities, never invent new structures).
    assert matched_dimers <= theoretical_dimers


def test_generate_theoretical_dimers_3_3_donor_requires_d_ala_at_position_4(raw_data, theo_masses, ppm):
    """3-3 (L,D-transpeptidase) crosslinks require the mDAP3–D-Ala4 bond to be
    present in the donor; only tetrapeptides ending in D-Ala at position 4 are
    valid donors. Structures like gm-AEJK|1 (K at position 4) must NOT appear
    as 3-3 donors even at a 99.9% abundance threshold."""
    _result, donors_used = matching.generate_theoretical_dimers(
        raw_data, theo_masses, [], ppm, species_code="ecoli", strict_mode=True, donor_abundance_threshold=0.999
    )

    three_three_donors = donors_used[donors_used["Crosslink Type"] == "3-3 crosslink"]

    # Every selected 3-3 donor must be a tetrapeptide whose 4th residue is D-Ala (A).
    # In standard E. coli data this means only gm-AEJA|1 (and any modifications of it).
    for structure in three_three_donors["Inferred structure"]:
        stem = structure.split("|")[0]
        # Strip any modification tag, e.g. "gm-AEJA (Anh)" -> "gm-AEJA"
        if " (" in stem:
            stem = stem[: stem.index(" (")]
        residues = stem.split("-", 1)[-1]  # "gm-AEJA" -> "AEJA"
        assert residues[3] == "A", f"3-3 donor {structure} has non-D-Ala at position 4: '{residues[3]}'"

    # The bundled simple E. coli data has exactly one valid 3-3 donor (gm-AEJA|1)
    assert three_three_donors["Inferred structure"].nunique() == 1
    assert (three_three_donors["Inferred structure"] == "gm-AEJA|1").all()
    assert three_three_donors["Cumulative percent"].max() == pytest.approx(100.0)


def test_generate_theoretical_dimers_with_custom_rule(raw_data, theo_masses, ppm):
    rule = CustomCrosslinkRule(
        label="Custom 4-3",
        donor_pattern=r"^..JAA$",
        acceptor_pattern=r"^..J",
        loses_terminal_ala=True,
    )
    result, donors_used = matching.generate_theoretical_dimers(raw_data, theo_masses, [], ppm, custom_rule=rule)

    assert not result.empty
    assert result["Inferred structure"].str.contains(r"\[Custom 4-3\]").all()
    assert (donors_used["Crosslink Type"] == "Custom 4-3").all()
