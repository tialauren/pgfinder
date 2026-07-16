"""Biological validation tests using real experimental ftrs data.

These tests check that the dimer-matching code produces biologically correct
output for each species — correct crosslink types, correct donors, correct
structure notation — against real MS1 data rather than synthetic fixtures.

They are separate from the unit tests (which use constructed DataFrames) and
from the regression tests (which snapshot exact output). The goal here is to
catch biochemical rule errors: wrong crosslink types, donors with the wrong
residue at position 4, bridge constraints ignored, etc.

Files in tests/resources/data/:
  ftrs_cdiff.ftrs        C. difficile WT real MS1 data
  ftrs_bsubtilis.ftrs    B. subtilis (strain Rosa) real MS1 data
  ftrs_efaecalis.ftrs    E. faecalis WT real MS1 data

Minimal ground-truth mass databases in tests/resources/data/masses/:
  cdiff_biological_core.csv      g(-Ac)m-AEJ/AEJA/AEJAA only
  bsubtilis_biological_core.csv  gm-AEJ/AEJA/AEJAA + gm-AQK/AQKA/AQKAA
  efaecalis_biological_core.csv  gm-AQKAA, gm-AQKA, gm-AQK[D]/[D]A/[D]AA

Species not covered here:
  Fusobacterium  needs a Fusobacterium-specific mass database (E. coli masses
                 do not match the ftrs data within normal ppm tolerance).
  S. aureus      no ftrs data available yet.
"""

from pathlib import Path

import pytest

import pgfinder.pgio as pgio
from pgfinder import matching
from pgfinder.accurate_species_rules import parse_peptide_structure

DATA   = Path(__file__).parent / "resources" / "data"
MASSES = DATA / "masses"
PPM    = 10


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _crosslink_types_in(dimers):
    """Return the set of crosslink type tags found in dimer structure strings."""
    import re

    found = set()
    for s in dimers["Inferred structure"]:
        m = re.search(r"\[([^\]]+)\]\|2$", s)
        if m:
            found.add(m.group(1))
    return found


# ---------------------------------------------------------------------------
# C. difficile
# ---------------------------------------------------------------------------


class TestCDiff:
    @pytest.fixture(scope="class")
    def results(self):
        raw    = pgio.ftrs_reader(DATA / "ftrs_cdiff.ftrs")
        masses = pgio.theo_masses_reader(MASSES / "cdiff_biological_core.csv")
        dimers, donors = matching.generate_theoretical_dimers(
            raw, masses, [], PPM, species_code="cdiff"
        )
        return dimers, donors

    def test_both_crosslink_types_produced(self, results):
        dimers, _ = results
        types = _crosslink_types_in(dimers)
        assert "4-3" in types, "4-3 crosslinks missing from C. diff output"
        assert "3-3" in types, "3-3 crosslinks missing from C. diff output"

    def test_no_other_crosslink_types(self, results):
        dimers, _ = results
        assert _crosslink_types_in(dimers) == {"4-3", "3-3"}

    def test_correct_4_3_donor(self, results):
        """4-3 donor must be the C. diff pentapeptide (de-N-acetylated glycan)."""
        _, donors = results
        four_three = donors[donors["Crosslink Type"] == "4-3 crosslink"]
        assert not four_three.empty, "No 4-3 donor found"
        for structure in four_three["Inferred structure"]:
            info = parse_peptide_structure(structure)
            assert info is not None, f"Could not parse donor: {structure}"
            assert info.length == 5, f"4-3 donor must be pentapeptide, got length {info.length}: {structure}"
            assert info.has_mdap_pos3, f"4-3 donor must have mDAP at pos3: {structure}"
            assert info.sequence[3] == "A", f"4-3 donor must have D-Ala at pos4: {structure}"
            assert info.sequence[4] == "A", f"4-3 donor must have D-Ala at pos5: {structure}"

    def test_correct_3_3_donor(self, results):
        """3-3 donor must be the C. diff tetrapeptide with D-Ala at pos4."""
        _, donors = results
        three_three = donors[donors["Crosslink Type"] == "3-3 crosslink"]
        assert not three_three.empty, "No 3-3 donor found"
        for structure in three_three["Inferred structure"]:
            info = parse_peptide_structure(structure)
            assert info is not None
            assert info.length == 4, f"3-3 donor must be tetrapeptide, got {info.length}: {structure}"
            assert info.has_mdap_pos3, f"3-3 donor must have mDAP at pos3: {structure}"
            assert info.sequence[3] == "A", f"3-3 donor must have D-Ala at pos4 (LDT recognition): {structure}"

    def test_donor_trimming_in_dimer_notation(self, results):
        """4-3 dimer notation must show tetrapeptide (not pentapeptide) donor.
        3-3 dimer notation must show tripeptide (not tetrapeptide) donor."""
        dimers, _ = results
        for _, row in dimers.iterrows():
            structure = row["Inferred structure"]
            donor_part = structure.split("=")[0]
            info = parse_peptide_structure(donor_part + "|2")
            if "[4-3]" in structure:
                assert info.length == 4, f"4-3 donor in dimer must be trimmed to tetrapeptide: {structure}"
            elif "[3-3]" in structure:
                assert info.length == 3, f"3-3 donor in dimer must be trimmed to tripeptide: {structure}"

    def test_all_dimers_tagged_as_dimer(self, results):
        dimers, _ = results
        assert dimers["Inferred structure"].str.endswith("|2").all()

    def test_isobaric_pairs_present(self, results):
        """The 4-3 tetra-tri and 3-3 tri-tetra dimers are isobaric in C. diff
        (same stem peptide rules as E. coli). Both should appear."""
        dimers, _ = results
        structures = set(dimers["Inferred structure"])
        assert "g(-Ac)m-AEJA=g(-Ac)m-AEJ[4-3]|2" in structures
        assert "g(-Ac)m-AEJ=g(-Ac)m-AEJA[3-3]|2" in structures


# ---------------------------------------------------------------------------
# B. subtilis
# ---------------------------------------------------------------------------


class TestBSubtilis:
    """This sample (strain Rosa) has detectable tetrapeptide but no pentapeptide,
    so only 3-3 crosslinks can be validated from this dataset. That is
    biologically plausible (high D,D-carboxypeptidase activity trimming
    pentapeptides). The 4-3 rule is validated separately in unit tests."""

    @pytest.fixture(scope="class")
    def results(self):
        raw    = pgio.ftrs_reader(DATA / "ftrs_bsubtilis.ftrs")
        masses = pgio.theo_masses_reader(MASSES / "bsubtilis_biological_core.csv")
        dimers, donors = matching.generate_theoretical_dimers(
            raw, masses, [], PPM, species_code="bsubtilis"
        )
        return dimers, donors

    def test_3_3_crosslinks_produced(self, results):
        dimers, _ = results
        assert not dimers.empty, "No dimers produced for B. subtilis"
        assert "3-3" in _crosslink_types_in(dimers)

    def test_3_3_donor_has_d_ala_at_pos4(self, results):
        """LDT requires D-Ala at position 4 — no other residue is a valid substrate."""
        _, donors = results
        three_three = donors[donors["Crosslink Type"] == "3-3 crosslink"]
        assert not three_three.empty, "No 3-3 donor found in B. subtilis data"
        for structure in three_three["Inferred structure"]:
            info = parse_peptide_structure(structure)
            assert info is not None
            assert info.length == 4, f"3-3 donor must be tetrapeptide: {structure}"
            assert info.sequence[3] == "A", (
                f"3-3 donor must have D-Ala at pos4 — non-Ala residue would not be "
                f"cleaved by L,D-transpeptidase: {structure}"
            )

    def test_no_bridge_crosslinks_in_dap_type_sample(self, results):
        """B. subtilis DAP-type strains do not use pentaglycine or D-Asp bridges."""
        dimers, _ = results
        types = _crosslink_types_in(dimers)
        assert "4-3-bridge" not in types


# ---------------------------------------------------------------------------
# E. faecalis
# ---------------------------------------------------------------------------


class TestEFaecalis:
    @pytest.fixture(scope="class")
    def results(self):
        raw    = pgio.ftrs_reader(DATA / "ftrs_efaecalis.ftrs")
        masses = pgio.theo_masses_reader(MASSES / "efaecalis_biological_core.csv")
        dimers, donors = matching.generate_theoretical_dimers(
            raw, masses, [], PPM, species_code="efaecalis"
        )
        return dimers, donors

    def test_both_crosslink_types_produced(self, results):
        dimers, _ = results
        types = _crosslink_types_in(dimers)
        assert "4-3-bridge" in types, "4-3-bridge crosslinks missing from E. faecalis output"
        assert "3-3" in types, "3-3 crosslinks missing from E. faecalis output"

    def test_no_direct_4_3_crosslinks(self, results):
        """E. faecalis crosslinks only through the Ala₂ bridge — direct 4-3 must not appear."""
        dimers, _ = results
        assert "4-3" not in _crosslink_types_in(dimers), (
            "Bare 4-3 crosslink appeared in E. faecalis output — only 4-3-bridge is valid"
        )

    def test_4_3_bridge_donor_carries_ala2_bridge(self, results):
        """4-3-bridge donor must carry the di-L-Ala (Ala₂) bridge [AA]."""
        _, donors = results
        bridge_donors = donors[donors["Crosslink Type"] == "4-3 bridge crosslink"]
        assert not bridge_donors.empty, "No 4-3-bridge donor found"
        for structure in bridge_donors["Inferred structure"]:
            info = parse_peptide_structure(structure)
            assert info is not None
            assert info.length == 5, f"4-3-bridge donor must be pentapeptide: {structure}"
            assert info.has_lys_pos3, f"4-3-bridge donor must have Lys at pos3: {structure}"
            assert info.has_ala2_bridge, (
                f"4-3-bridge donor must carry the [AA] (di-L-Ala) bridge: {structure}"
            )
            assert info.sequence[3] == "A", f"D-Ala at pos4 required: {structure}"
            assert info.sequence[4] == "A", f"D-Ala at pos5 required: {structure}"

    def test_4_3_bridge_acceptors_carry_ala2_bridge(self, results):
        """Every 4-3-bridge dimer must reference an acceptor that has the [AA] (Ala₂) bridge.

        Dimer format: donor=acceptor[crosslink_type]|2
        e.g. gm-AQK[AA]A=gm-AQK[AA][4-3-bridge]|2
        Strip the trailing [crosslink_type]|2 to isolate the acceptor.
        """
        import re

        dimers, _ = results
        bridge_dimers = dimers[dimers["Inferred structure"].str.contains(r"\[4-3-bridge\]")]
        assert not bridge_dimers.empty
        for _, row in bridge_dimers.iterrows():
            structure = row["Inferred structure"]
            right = structure.split("=")[1]           # gm-AQK[AA][4-3-bridge]|2
            right = right.rsplit("|", 1)[0]           # gm-AQK[AA][4-3-bridge]
            acceptor_part = re.sub(r"\[[^\]]+\]$", "", right)  # gm-AQK[AA]
            assert "[AA]" in acceptor_part, (
                f"4-3-bridge acceptor must carry [AA] (di-L-Ala bridge): {structure}"
            )

    def test_3_3_donor_has_lys_pos3_and_d_ala_pos4(self, results):
        _, donors = results
        three_three = donors[donors["Crosslink Type"] == "3-3 crosslink"]
        assert not three_three.empty, "No 3-3 donor found"
        for structure in three_three["Inferred structure"]:
            info = parse_peptide_structure(structure)
            assert info is not None
            assert info.length == 4, f"3-3 donor must be tetrapeptide: {structure}"
            assert info.has_lys_pos3, f"3-3 donor must have Lys at pos3: {structure}"
            assert info.sequence[3] == "A", f"3-3 donor must have D-Ala at pos4: {structure}"
            assert not info.has_ala2_bridge, (
                f"3-3 donor must not have [AA] bridge — Lys3 epsilon-NH2 already occupied: {structure}"
            )

    def test_all_dimers_tagged_as_dimer(self, results):
        dimers, _ = results
        assert dimers["Inferred structure"].str.endswith("|2").all()
