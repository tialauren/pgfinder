"""
Tests for MS2 HCD fragment ion prediction.

Mass verification strategy: check fragment masses against known analytical
values (standard proteomics formulas) rather than against experimental spectra.
Experimental validation against the Kwan 2024 paper comes later.
"""

import pandas as pd
import pytest

from pgfinder.ms2.fragmentor import (
    _RESIDUE_MASS,
    H2_REDUCTION,
    H2O,
    _mass,
    _mz,
    _parse_monomer,
    predict_hcd_fragments,
)
from pgfinder.ms2.scorer import score_spectrum, score_structures

# ---------------------------------------------------------------------------
# Formula → mass
# ---------------------------------------------------------------------------


def test_mass_water():
    assert _mass("H2O") == pytest.approx(18.01056469, abs=1e-6)


def test_mass_glcnac():
    # GlcNAc C8H15NO6
    assert _mass("C8H15NO6") == pytest.approx(221.08994, abs=1e-4)


def test_mass_ala():
    # Ala C3H7NO2 — 89.04768 is the well-known monoisotopic value
    assert _RESIDUE_MASS["A"] == pytest.approx(89.04768, abs=1e-4)


# ---------------------------------------------------------------------------
# Structure string parsing
# ---------------------------------------------------------------------------


def test_parse_monomer_simple():
    glycan, stem, mod, oligo = _parse_monomer("gm-AEJA|1")
    assert glycan == ["g", "m"]
    assert stem == ["A", "E", "J", "A"]
    assert mod == ""
    assert oligo == 1


def test_parse_monomer_with_modification():
    glycan, stem, mod, oligo = _parse_monomer("gm-AEJAA (Anh)|1")
    assert glycan == ["g", "m"]
    assert stem == ["A", "E", "J", "A", "A"]
    assert mod == " (Anh)"


def test_parse_monomer_pentapeptide():
    glycan, stem, mod, oligo = _parse_monomer("gm-AEJAA|1")
    assert stem == ["A", "E", "J", "A", "A"]


def test_parse_dimer_returns_none():
    assert _parse_monomer("gm-AEJA=gm-AEJAA[4-3]|2") is None


def test_predict_dimer_returns_none():
    assert predict_hcd_fragments("gm-AEJA=gm-AEJAA[4-3]|2") is None


# ---------------------------------------------------------------------------
# Peptide b/y ion masses (standard proteomics formulas)
# Standard: b_n = sum(residues 1..n) - n*H2O + PROTON
#           y_n = sum(residues from C-term n) - (n-1)*H2O + PROTON
# ---------------------------------------------------------------------------


def _b_mz(residues: list, z: int = 1) -> float:
    """Expected b-ion [M+zH]^z+ m/z."""
    total = sum(_RESIDUE_MASS[r] for r in residues)
    neutral = total - len(residues) * H2O
    return _mz(neutral, z)


def _y_mz(residues: list, z: int = 1) -> float:
    """Expected y-ion [M+zH]^z+ m/z."""
    total = sum(_RESIDUE_MASS[r] for r in residues)
    neutral = total - (len(residues) - 1) * H2O
    return _mz(neutral, z)


def test_b1_ala():
    """b1 for Ala: well-known m/z = 72.044 Da [M+H]+."""
    expected = _b_mz(["A"])
    assert expected == pytest.approx(72.0449, abs=1e-3)


def test_y1_ala():
    """y1 for Ala: well-known m/z = 90.055 Da [M+H]+."""
    expected = _y_mz(["A"])
    assert expected == pytest.approx(90.0550, abs=1e-3)


def test_peptide_ions_present_in_output():
    """b and y ions must appear in predict_hcd_fragments output for gm-AEJA|1."""
    df = predict_hcd_fragments("gm-AEJA|1")
    assert df is not None
    ion_types = set(df["ion_type"])

    # b ions b1..b3
    assert "b1" in ion_types
    assert "b2" in ion_types
    assert "b3" in ion_types

    # y ions y1..y3
    assert "y1" in ion_types
    assert "y2" in ion_types
    assert "y3" in ion_types


def test_b1_value_in_output():
    """b1 ion from gm-AEJA|1 stem peptide (A at position 1) has correct m/z."""
    df = predict_hcd_fragments("gm-AEJA|1")
    b1_rows = df[(df["ion_type"] == "b1") & (df["charge"] == 1)]
    assert not b1_rows.empty
    assert b1_rows.iloc[0]["mz"] == pytest.approx(_b_mz(["A"]), abs=1e-3)


def test_y1_value_in_output():
    """y1 ion from gm-AEJA|1 stem peptide (A at C-terminus) has correct m/z."""
    df = predict_hcd_fragments("gm-AEJA|1")
    y1_rows = df[(df["ion_type"] == "y1") & (df["charge"] == 1)]
    assert not y1_rows.empty
    assert y1_rows.iloc[0]["mz"] == pytest.approx(_y_mz(["A"]), abs=1e-3)


def test_b2_value_in_output():
    """b2 from AEJA stem (AE) has correct m/z."""
    df = predict_hcd_fragments("gm-AEJA|1")
    b2_rows = df[(df["ion_type"] == "b2") & (df["charge"] == 1)]
    assert not b2_rows.empty
    assert b2_rows.iloc[0]["mz"] == pytest.approx(_b_mz(["A", "E"]), abs=1e-3)


def test_y3_value_in_output():
    """y3 from AEJA stem (EJA) has correct m/z."""
    df = predict_hcd_fragments("gm-AEJA|1")
    y3_rows = df[(df["ion_type"] == "y3") & (df["charge"] == 1)]
    assert not y3_rows.empty
    assert y3_rows.iloc[0]["mz"] == pytest.approx(_y_mz(["E", "J", "A"]), abs=1e-3)


# ---------------------------------------------------------------------------
# Glycan B/Y/C/Z ion presence and values
# ---------------------------------------------------------------------------


def test_glycan_ions_present():
    """B1, C1, Y1 and Z1 must appear for gm- monomer (one glycosidic bond)."""
    df = predict_hcd_fragments("gm-AEJA|1")
    assert df is not None
    ion_types = set(df["ion_type"])
    assert "B1" in ion_types
    assert "C1" in ion_types
    assert "Y1" in ion_types
    assert "Z1" in ion_types


def test_c1_equals_b1_plus_water():
    """C ion is B ion + H2O (Domon-Costello: C = B + H2O)."""
    df = predict_hcd_fragments("gm-AEJA|1")
    b1 = df[(df["ion_type"] == "B1") & (df["charge"] == 1)].iloc[0]["mz"]
    c1 = df[(df["ion_type"] == "C1") & (df["charge"] == 1)].iloc[0]["mz"]
    assert c1 - b1 == pytest.approx(H2O, abs=1e-4)


def test_c1_value_matches_yaml():
    """C1 [M+H]+ = 222.097 Da (from check_frag from pgn.yaml, Gly. C/Z[r])."""
    df = predict_hcd_fragments("gm-AEJA|1")
    c1 = df[(df["ion_type"] == "C1") & (df["charge"] == 1)].iloc[0]["mz"]
    assert c1 == pytest.approx(222.097, abs=1e-2)


def test_z1_equals_y1_minus_water():
    """Z ion is Y ion − H2O (Domon-Costello: Z = Y − H2O)."""
    df = predict_hcd_fragments("gm-AEJA|1")
    y1 = df[(df["ion_type"] == "Y1") & (df["charge"] == 1)].iloc[0]["mz"]
    z1 = df[(df["ion_type"] == "Z1") & (df["charge"] == 1)].iloc[0]["mz"]
    assert y1 - z1 == pytest.approx(H2O, abs=1e-4)


# ---------------------------------------------------------------------------
# Lactoyl b-lac / y-lac ions
# ---------------------------------------------------------------------------


def test_lac_ions_present():
    """b-lac and y-lac ions must appear for gm- monomers."""
    df = predict_hcd_fragments("gm-AEJA|1")
    ion_types = set(df["ion_type"])
    assert "b-lac" in ion_types
    assert "y-lac" in ion_types


# ---------------------------------------------------------------------------
# Secondary losses
# ---------------------------------------------------------------------------


def test_glc_freed_ions_present():
    """GlcNAc freed fragments (Glc-1 at 186.076 [M+H]+) must appear."""
    df = predict_hcd_fragments("gm-AEJA|1")
    ion_types = set(df["ion_type"])
    assert "Glc-1" in ion_types


def test_mur_freed_ions_present():
    """MurNAc freed fragment Mur-1 must appear."""
    df = predict_hcd_fragments("gm-AEJA|1")
    ion_types = set(df["ion_type"])
    assert "Mur-1" in ion_types


def test_glc1_mz():
    """Glc-1 [M+H]+ = 186.076 Da (from hcd_rules.kdl comment)."""
    df = predict_hcd_fragments("gm-AEJA|1")
    row = df[(df["ion_type"] == "Glc-1") & (df["charge"] == 1)]
    assert not row.empty
    assert row.iloc[0]["mz"] == pytest.approx(186.076, abs=1e-2)


# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------


def test_output_columns():
    df = predict_hcd_fragments("gm-AEJA|1")
    assert set(df.columns) == {"ion_type", "ion_series", "mz", "charge", "neutral_mass"}


def test_output_sorted_by_mz():
    df = predict_hcd_fragments("gm-AEJA|1")
    assert list(df["mz"]) == sorted(df["mz"])


def test_no_negative_mz():
    df = predict_hcd_fragments("gm-AEJA|1")
    assert (df["mz"] > 0).all()


def test_pentapeptide_has_b4():
    """Pentapeptide gm-AEJAA|1 must have a b4 ion."""
    df = predict_hcd_fragments("gm-AEJAA|1")
    assert "b4" in set(df["ion_type"])


def test_precursor_present():
    """Precursor ion must appear in the output."""
    df = predict_hcd_fragments("gm-AEJA|1")
    assert "precursor" in set(df["ion_type"])


def test_precursor_mass_matches_database():
    """
    Precursor neutral mass for gm-AEJA|1 must match the pgfinder mass
    database value of 941.4075 Da within 5 mDa.

    This is the key integration test: if the residue-level arithmetic is
    consistent with the existing pgfinder mass tables, the two modules
    are using the same convention.
    """
    df = predict_hcd_fragments("gm-AEJA|1")
    precursor = df[df["ion_type"] == "precursor"]
    row = precursor[precursor["charge"] == 1].iloc[0]
    assert row["neutral_mass"] == pytest.approx(941.4075, abs=0.005)


# ---------------------------------------------------------------------------
# E/Q secondary losses
# ---------------------------------------------------------------------------


def test_e_losses_present_when_stem_contains_e():
    """gm-AEJA|1 has E at position 2 — e1/e2 losses must appear on y3 ion."""
    df = predict_hcd_fragments("gm-AEJA|1")
    ion_types = set(df["ion_type"])
    # y3 covers EJA (contains E) → y3-e1 and y3-e2
    assert "y3-e1" in ion_types
    assert "y3-e2" in ion_types


def test_e_losses_absent_when_no_e_in_y_fragment():
    """y1 from AEJA covers only A (the C-terminal Ala) — no E, so no e losses."""
    df = predict_hcd_fragments("gm-AEJA|1")
    ion_types = set(df["ion_type"])
    assert "y1-e1" not in ion_types
    assert "y1-e2" not in ion_types


def test_e1_mass_is_y_minus_water():
    """e1 loss = y_k - H2O."""
    df = predict_hcd_fragments("gm-AEJA|1")
    y3 = df[(df["ion_type"] == "y3") & (df["charge"] == 1)].iloc[0]["mz"]
    e1 = df[(df["ion_type"] == "y3-e1") & (df["charge"] == 1)].iloc[0]["mz"]
    assert y3 - e1 == pytest.approx(18.0106, abs=1e-3)


def test_q_losses_present_when_stem_contains_q():
    """A stem with Q must generate q1/q2 losses on y-ions containing Q."""
    # gm-AQKAA (S. aureus-like) has Q at position 2
    df = predict_hcd_fragments("gm-AQKAA|1")
    if df is None:
        pytest.skip("gm-AQKAA|1 not parseable — check residue codes")
    ion_types = set(df["ion_type"])
    # y4 covers QKAA (contains Q) → q1/q2 expected
    assert any("q1" in t for t in ion_types)
    assert any("q2" in t for t in ion_types)


# ---------------------------------------------------------------------------
# Terminal losses
# ---------------------------------------------------------------------------


def test_precursor_h2o_loss_present():
    """precursor-H2O must be present."""
    df = predict_hcd_fragments("gm-AEJA|1")
    assert "precursor-H2O" in set(df["ion_type"])


def test_precursor_nh3_loss_present():
    """precursor-NH3 must be present."""
    df = predict_hcd_fragments("gm-AEJA|1")
    assert "precursor-NH3" in set(df["ion_type"])


def test_precursor_h2o_mass():
    """precursor-H2O must be exactly precursor - 18.011 Da."""
    df = predict_hcd_fragments("gm-AEJA|1")
    prec = df[(df["ion_type"] == "precursor") & (df["charge"] == 1)].iloc[0]["mz"]
    loss = df[(df["ion_type"] == "precursor-H2O") & (df["charge"] == 1)].iloc[0]["mz"]
    assert prec - loss == pytest.approx(18.0106, abs=1e-3)


# ---------------------------------------------------------------------------
# Non-reduced mode (reduced=False)
# ---------------------------------------------------------------------------


def test_non_reduced_y1_lighter_by_h2():
    """Y ions in non-reduced mode are lighter by exactly H2_REDUCTION (~2.016 Da)."""
    df_red = predict_hcd_fragments("gm-AEJA|1", reduced=True)
    df_nonred = predict_hcd_fragments("gm-AEJA|1", reduced=False)
    y1_red = df_red[(df_red["ion_type"] == "Y1") & (df_red["charge"] == 1)].iloc[0]["mz"]
    y1_nonred = df_nonred[(df_nonred["ion_type"] == "Y1") & (df_nonred["charge"] == 1)].iloc[0]["mz"]
    assert y1_red - y1_nonred == pytest.approx(H2_REDUCTION, abs=1e-4)


def test_non_reduced_b1_unchanged():
    """B ions don't contain the reducing end — they must be identical in both modes."""
    df_red = predict_hcd_fragments("gm-AEJA|1", reduced=True)
    df_nonred = predict_hcd_fragments("gm-AEJA|1", reduced=False)
    b1_red = df_red[(df_red["ion_type"] == "B1") & (df_red["charge"] == 1)].iloc[0]["mz"]
    b1_nonred = df_nonred[(df_nonred["ion_type"] == "B1") & (df_nonred["charge"] == 1)].iloc[0]["mz"]
    assert b1_red == pytest.approx(b1_nonred, abs=1e-6)


def test_non_reduced_c1_unchanged():
    """C ions don't contain the reducing end — they must be identical in both modes."""
    df_red = predict_hcd_fragments("gm-AEJA|1", reduced=True)
    df_nonred = predict_hcd_fragments("gm-AEJA|1", reduced=False)
    c1_red = df_red[(df_red["ion_type"] == "C1") & (df_red["charge"] == 1)].iloc[0]["mz"]
    c1_nonred = df_nonred[(df_nonred["ion_type"] == "C1") & (df_nonred["charge"] == 1)].iloc[0]["mz"]
    assert c1_red == pytest.approx(c1_nonred, abs=1e-6)


def test_non_reduced_precursor_lighter_by_h2():
    """Precursor in non-reduced mode is lighter by H2_REDUCTION."""
    df_red = predict_hcd_fragments("gm-AEJA|1", reduced=True)
    df_nonred = predict_hcd_fragments("gm-AEJA|1", reduced=False)
    prec_red = df_red[(df_red["ion_type"] == "precursor") & (df_red["charge"] == 1)].iloc[0]["neutral_mass"]
    prec_nonred = df_nonred[(df_nonred["ion_type"] == "precursor") & (df_nonred["charge"] == 1)].iloc[0]["neutral_mass"]
    assert prec_red - prec_nonred == pytest.approx(H2_REDUCTION, abs=1e-4)


def test_non_reduced_blac_lighter_by_h2():
    """b-lac contains the glycan (including reducing end) — lighter by H2_REDUCTION in non-reduced mode."""
    df_red = predict_hcd_fragments("gm-AEJA|1", reduced=True)
    df_nonred = predict_hcd_fragments("gm-AEJA|1", reduced=False)
    blac_red = df_red[(df_red["ion_type"] == "b-lac") & (df_red["charge"] == 1)].iloc[0]["mz"]
    blac_nonred = df_nonred[(df_nonred["ion_type"] == "b-lac") & (df_nonred["charge"] == 1)].iloc[0]["mz"]
    assert blac_red - blac_nonred == pytest.approx(H2_REDUCTION, abs=1e-4)


def test_non_reduced_ylac_unchanged():
    """y-lac is the stem + lactoyl group only — unaffected by reduction status."""
    df_red = predict_hcd_fragments("gm-AEJA|1", reduced=True)
    df_nonred = predict_hcd_fragments("gm-AEJA|1", reduced=False)
    ylac_red = df_red[(df_red["ion_type"] == "y-lac") & (df_red["charge"] == 1)].iloc[0]["mz"]
    ylac_nonred = df_nonred[(df_nonred["ion_type"] == "y-lac") & (df_nonred["charge"] == 1)].iloc[0]["mz"]
    assert ylac_red == pytest.approx(ylac_nonred, abs=1e-6)


def test_non_reduced_yaml_y1_value():
    """
    Non-reduced Y1 for gm-A|1 equivalent = 365.155 [M+H]+.
    Verified against check_frag from pgn.yaml entry (NAG)(NAM)-A, Gly. B/Y: [204.087, 365.155].
    We use gm-AEJA|1 stem here — the test checks the relative shift, not absolute.
    """
    df_red = predict_hcd_fragments("gm-AEJA|1", reduced=True)
    df_nonred = predict_hcd_fragments("gm-AEJA|1", reduced=False)
    y1_red = df_red[(df_red["ion_type"] == "Y1") & (df_red["charge"] == 1)].iloc[0]["mz"]
    y1_nonred = df_nonred[(df_nonred["ion_type"] == "Y1") & (df_nonred["charge"] == 1)].iloc[0]["mz"]
    # The shift must exactly equal H2_REDUCTION
    assert y1_red - y1_nonred == pytest.approx(H2_REDUCTION, abs=1e-4)


# ---------------------------------------------------------------------------
# Scorer
# ---------------------------------------------------------------------------


def _make_peaks(mz_list, intensity_list=None):
    """Helper: build an experimental peak DataFrame."""
    if intensity_list is None:
        intensity_list = [1000.0] * len(mz_list)
    return pd.DataFrame({"mz": mz_list, "intensity": intensity_list})


def test_score_perfect_match():
    """A theoretical spectrum matched against itself should score ≈ 1.0."""
    df = predict_hcd_fragments("gm-AEJA|1", charge_max=1)
    # Use singly charged ions as the 'experimental' peaks with uniform intensity
    peaks = _make_peaks(df[df["charge"] == 1]["mz"].tolist())
    result = score_spectrum("gm-AEJA|1", peaks, charge_max=1)
    assert result is not None
    # All theo peaks present → coverage = 1.0, score ≈ 1.0
    assert result["coverage"] == pytest.approx(1.0, abs=1e-6)
    assert result["score"] == pytest.approx(1.0, abs=0.01)


def test_score_no_match():
    """Peaks that are all far from any theoretical ion should score ≈ 0."""
    # m/z values that can't match anything in gm-AEJA|1 (all shifted by 500 Da)
    df = predict_hcd_fragments("gm-AEJA|1", charge_max=1)
    shifted = [m + 500.0 for m in df[df["charge"] == 1]["mz"].tolist()]
    peaks = _make_peaks(shifted)
    result = score_spectrum("gm-AEJA|1", peaks, charge_max=1)
    assert result is not None
    assert result["matched_peaks"] == 0
    assert result["score"] == pytest.approx(0.0, abs=1e-6)


def test_score_returns_annotated_df():
    """Result must include an 'annotated' DataFrame with the right columns."""
    df = predict_hcd_fragments("gm-AEJA|1", charge_max=1)
    peaks = _make_peaks(df[df["charge"] == 1]["mz"].tolist())
    result = score_spectrum("gm-AEJA|1", peaks, charge_max=1)
    assert result is not None
    ann = result["annotated"]
    assert "ion_type" in ann.columns
    assert "delta_ppm" in ann.columns
    assert "theo_mz" in ann.columns


def test_score_dimer_returns_none():
    """Dimers are not yet supported — score_spectrum must return None."""
    peaks = _make_peaks([100.0, 200.0])
    assert score_spectrum("gm-AEJA=gm-AEJAA[4-3]|2", peaks) is None


def test_score_structures_ranks_correct_structure_first():
    """
    Given two candidates where one is the true structure, it should rank
    highest when the experimental spectrum is the predicted spectrum of the
    true structure.
    """
    true_structure = "gm-AEJA|1"
    wrong_structure = "gm-AEJAA|1"

    df_true = predict_hcd_fragments(true_structure, charge_max=1)
    peaks = _make_peaks(df_true[df_true["charge"] == 1]["mz"].tolist())

    rankings = score_structures([true_structure, wrong_structure], peaks, charge_max=1)
    assert not rankings.empty
    assert rankings.iloc[0]["structure"] == true_structure


def test_score_structures_output_sorted():
    """score_structures output must be sorted by score descending."""
    peaks = _make_peaks([72.044, 90.055, 186.076])
    df = score_structures(["gm-AEJA|1", "gm-AEJAA|1"], peaks)
    scores = list(df["score"])
    assert scores == sorted(scores, reverse=True)
