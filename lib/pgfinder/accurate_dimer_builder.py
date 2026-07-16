"""
Accurate dimer builder using biochemically correct species rules.

Finds donors from detected monomers (MS1 results).
Finds acceptors from theoretical database (all possibilities).
"""

import logging
from decimal import Decimal
from typing import Dict, List, Optional, Tuple

import pandas as pd

from pgfinder import COLUMNS
from pgfinder.accurate_species_rules import (
    AccurateSpeciesRules,
    CustomCrosslinkRule,
    match_custom_acceptor,
    match_custom_donor,
    parse_peptide_structure,
    split_modification_tag,
)
from pgfinder.logs.logs import LOGGER_NAME

LOGGER = logging.getLogger(LOGGER_NAME)

# Use pgfinder columns
COLUMNS = COLUMNS["pgfinder"]

# Display labels for donors_used output. Built-in crosslink type codes like
# "4-3" are misread as dates (April 3) by Excel when a CSV is double-clicked —
# appending " crosslink" makes the value unambiguously non-date. Custom rule
# labels (e.g. "Custom 4-3") start with a word so they pass through unchanged.
_CROSSLINK_TYPE_LABELS: dict = {
    "4-3": "4-3 crosslink",
    "3-3": "3-3 crosslink",
    "4-3-bridge": "4-3 bridge crosslink",
}


class AccurateDimerBuilder:
    """
    Dimer builder that uses biochemically accurate crosslinking rules.

    Donors: Selected from detected monomers (MS1 results)
    Acceptors: Selected from theoretical database (all structures)
    """

    def __init__(
        self,
        species_code: Optional[str] = None,
        custom_rule: Optional[CustomCrosslinkRule] = None,
        strict_mode: bool = True,
        donor_abundance_threshold: float = 0.9,
    ):
        """
        Initialize with species-specific (or user-supplied custom) rules.

        Parameters
        ----------
        species_code : str, optional
            Species identifier (e.g., 'ecoli', 'saureus', 'efaecalis'). Exactly
            one of species_code/custom_rule must be given.
        custom_rule : CustomCrosslinkRule, optional
            User-defined donor/acceptor regex rule for a species not covered
            by the built-in species list. Exactly one of species_code/custom_rule
            must be given.
        strict_mode : bool
            If True: Enforce exact literature rules (default). Ignored when
            custom_rule is given.
            If False: Allow bridge variants for novel discovery
        donor_abundance_threshold : float
            Fraction (0.0-1.0) of the *structurally eligible* donor pool's
            cumulative intensity to cover when selecting donors, per crosslink
            type - see `_find_eligible_donors_within_threshold`. Defaults to
            0.9 (90%). A donor pool's most abundant donor is always included
            regardless of how low this is set.
        """
        if (species_code is None) == (custom_rule is None):
            raise ValueError("Exactly one of species_code or custom_rule must be given")

        self.custom_rule = custom_rule
        self.strict_mode = strict_mode
        self.donor_abundance_threshold = donor_abundance_threshold

        if custom_rule is not None:
            self.species_code = None
            self.species_info = None
            self.crosslink_types = [custom_rule.label]
            LOGGER.info(f"AccurateDimerBuilder initialized for custom rule '{custom_rule.label}'")
        else:
            self.species_code = species_code.lower()
            self.species_info = AccurateSpeciesRules.get_species_info(self.species_code)
            self.crosslink_types = AccurateSpeciesRules.get_crosslink_types(self.species_code)

            mode_str = "STRICT (literature rules)" if strict_mode else "PERMISSIVE (novel discovery)"
            LOGGER.info(f"AccurateDimerBuilder initialized for {self.species_info['name']}")
            LOGGER.info(f"Mode: {mode_str}")
        LOGGER.info(f"Supported crosslink types: {', '.join(self.crosslink_types)}")

    def _is_valid_donor(self, structure: str, crosslink_type: str) -> bool:
        if self.custom_rule is not None:
            return match_custom_donor(structure, self.custom_rule)
        return AccurateSpeciesRules.is_valid_donor(structure, self.species_code, crosslink_type)

    def _is_valid_acceptor(self, structure: str, crosslink_type: str) -> bool:
        if self.custom_rule is not None:
            return match_custom_acceptor(structure, self.custom_rule)
        return AccurateSpeciesRules.is_valid_acceptor(structure, self.species_code, crosslink_type, self.strict_mode)

    def _get_lost_residues(self, crosslink_type: str) -> List[str]:
        if self.custom_rule is not None:
            return ["A"] if self.custom_rule.loses_terminal_ala else []
        return AccurateSpeciesRules.get_lost_residues(crosslink_type)

    def build_dimers_from_detected_monomers(
        self, matched_monomers_df: pd.DataFrame, theoretical_database_df: pd.DataFrame, columns: dict = COLUMNS
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Generate theoretical dimers.

        Donors: From detected monomers (MS1 results), selected by cumulative
            abundance - see `_find_eligible_donors_within_threshold`.
        Acceptors: From theoretical database - all possible structures

        Parameters
        ----------
        matched_monomers_df : pd.DataFrame
            DataFrame containing matched monomers with abundances (MS1 results)
        theoretical_database_df : pd.DataFrame
            DataFrame containing all theoretical structures from database
        columns : dict
            Column name mapping

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            (theoretical dimer structures and masses, donors actually used to
            build them - one row per donor per crosslink type, with its
            intensity and share of that crosslink type's eligible donor pool)
        """
        empty_dimers = pd.DataFrame(columns=[columns["inferred"]["structure"], columns["inferred"]["mass"]])
        empty_donors = pd.DataFrame(
            columns=[
                "Crosslink Type",
                columns["inferred"]["structure"],
                columns["input"]["intensity"],
                "Percent of eligible total",
                "Cumulative percent",
            ]
        )

        # Filter to only monomers (oligomerization state = 1).
        # Accept both "gm-AEJAA|1" (standard pgfinder format) and "gm-AEJAA"
        # (no suffix — some third-party mass databases omit the |N notation).
        # Explicitly exclude |2 dimers in either case.
        def _is_monomer(s):
            if not isinstance(s, str):
                return False
            return s.endswith("|1") or ("|" not in s and len(s) > 0)

        monomers_only = matched_monomers_df[
            matched_monomers_df[columns["inferred"]["structure"]].apply(_is_monomer)
        ].copy()

        if monomers_only.empty:
            LOGGER.warning("No monomers found in matched data. Cannot build dimers.")
            return empty_dimers, empty_donors

        LOGGER.info(f"Building dimers from {len(monomers_only)} detected monomers")

        all_dimers = []
        donors_used = []

        # Process each crosslink type
        for crosslink_type in self.crosslink_types:
            LOGGER.info(f"Processing {crosslink_type} crosslinks")

            # Find the eligible donors covering donor_abundance_threshold of the
            # eligible pool's abundance, from detected monomers
            donors = self._find_eligible_donors_within_threshold(monomers_only, crosslink_type, columns)

            if not donors:
                LOGGER.warning(f"No valid {crosslink_type} donor found in detected monomers")
                continue

            LOGGER.info(
                f"Selected {len(donors)} {crosslink_type} donor(s) covering "
                f"{donors[-1]['cumulative_percent']:.1f}% of the eligible pool"
            )

            # Find all valid ACCEPTORS from theoretical database (shared across all donors)
            acceptors = self._find_valid_acceptors_from_database(theoretical_database_df, crosslink_type, columns)

            if acceptors.empty:
                LOGGER.warning(f"No valid {crosslink_type} acceptors found in database")
                continue

            LOGGER.info(f"Found {len(acceptors)} valid {crosslink_type} acceptors in database")
            for i, (_, acceptor) in enumerate(acceptors.iterrows(), 1):
                LOGGER.debug(f"  Acceptor {i}: {acceptor['structure']}")

            for donor in donors:
                LOGGER.info(
                    f"  Donor: {donor['structure']} (intensity: {donor['intensity']:.2e}, "
                    f"{donor['percent_of_eligible_total']:.1f}% of eligible pool)"
                )

                # Generate dimers for this donor
                dimers = self._generate_dimers(
                    donor_structure=donor["structure"],
                    donor_mass=donor["mass"],
                    acceptors_df=acceptors,
                    crosslink_type=crosslink_type,
                    columns=columns,
                )
                all_dimers.extend(dimers)

                donors_used.append(
                    {
                        "Crosslink Type": _CROSSLINK_TYPE_LABELS.get(crosslink_type, crosslink_type),
                        columns["inferred"]["structure"]: donor["structure"],
                        columns["input"]["intensity"]: donor["intensity"],
                        "Percent of eligible total": donor["percent_of_eligible_total"],
                        "Cumulative percent": donor["cumulative_percent"],
                    }
                )

        dimers_df = pd.DataFrame(all_dimers) if all_dimers else empty_dimers
        donors_used_df = pd.DataFrame(donors_used) if donors_used else empty_donors

        if all_dimers:
            LOGGER.info(f"Generated {len(all_dimers)} theoretical dimers in total")
        else:
            LOGGER.warning("No dimers generated for any crosslink type")

        return dimers_df, donors_used_df

    def _find_eligible_donors_within_threshold(
        self, monomers_df: pd.DataFrame, crosslink_type: str, columns: dict
    ) -> List[dict]:
        """
        Find the most abundant structurally-eligible donors, covering
        `self.donor_abundance_threshold` of the *eligible donor pool's*
        cumulative intensity (not all detected monomers - an abundant but
        structurally ineligible monomer must not shrink the budget available
        to actual donor candidates).

        Always includes at least the single most abundant eligible donor,
        regardless of how low the threshold is; a threshold of 1.0 includes
        every eligible donor.

        The same structure can appear as multiple rows in `monomers_df` (e.g.
        matched at different retention times or charge states) - these are
        the same underlying donor candidate, so their intensities are summed
        before ranking, rather than treating each row as a separate donor
        (which would otherwise generate duplicate dimers).

        Searches in DETECTED monomers from MS1.

        Returns
        -------
        List[dict]
            One dict per selected donor, each with keys: structure, mass
            (Decimal), intensity, percent_of_eligible_total, cumulative_percent
            (the latter two as percentages, e.g. 37.5 not 0.375). Empty if no
            monomer is structurally eligible.
        """
        structure_column = columns["inferred"]["structure"]
        eligible = monomers_df[monomers_df[structure_column].apply(lambda s: self._is_valid_donor(s, crosslink_type))]
        if eligible.empty:
            return []

        eligible = eligible.groupby(structure_column, as_index=False).agg(
            **{
                columns["input"]["intensity"]: (columns["input"]["intensity"], "sum"),
                columns["inferred"]["mass"]: (columns["inferred"]["mass"], "first"),
            }
        )
        eligible = eligible.sort_values(by=columns["input"]["intensity"], ascending=False)
        eligible_total = eligible[columns["input"]["intensity"]].sum()

        donors = []
        cumulative = 0.0
        for _, row in eligible.iterrows():
            structure = row[columns["inferred"]["structure"]]
            intensity = row[columns["input"]["intensity"]]
            cumulative += intensity

            donors.append(
                {
                    "structure": structure,
                    "mass": Decimal(str(row[columns["inferred"]["mass"]])),
                    "intensity": intensity,
                    "percent_of_eligible_total": (intensity / eligible_total) * 100,
                    "cumulative_percent": (cumulative / eligible_total) * 100,
                }
            )

            info = parse_peptide_structure(structure)
            if info:
                LOGGER.debug(
                    f"  Valid donor found: {structure} (sequence: {''.join(info.sequence)}, length: {info.length})"
                )

            if cumulative / eligible_total >= self.donor_abundance_threshold:
                break

        return donors

    def _find_valid_acceptors_from_database(
        self, database_df: pd.DataFrame, crosslink_type: str, columns: dict
    ) -> pd.DataFrame:
        """
        Find all valid acceptors from the THEORETICAL DATABASE.

        This is different from finding donors - acceptors come from ALL
        possible structures in the database, not just detected ones.
        """
        acceptors = []
        structure_column = columns["inferred"]["structure"]
        mass_column = columns["inferred"]["mass"]

        for _, row in database_df.iterrows():
            structure = row[structure_column]

            # Check if this structure is a valid acceptor
            if self._is_valid_acceptor(structure, crosslink_type):
                acceptors.append({"structure": structure, "mass": row[mass_column]})

                # Log acceptor details
                info = parse_peptide_structure(structure)
                if info:
                    LOGGER.debug(f"  Valid acceptor: {structure}")
                    if self.custom_rule is None and crosslink_type == "4-3-bridge" and self.species_code == "saureus":
                        LOGGER.debug(f"    Bridge length: {info.glycine_bridge_length}")
                    elif (
                        self.custom_rule is None and crosslink_type == "4-3-bridge" and self.species_code == "efaecalis"
                    ):
                        LOGGER.debug(f"    Has Ala2 bridge: {info.has_ala2_bridge}")

        if not acceptors:
            return pd.DataFrame()

        return pd.DataFrame(acceptors)

    def _generate_dimers(
        self, donor_structure: str, donor_mass: Decimal, acceptors_df: pd.DataFrame, crosslink_type: str, columns: dict
    ) -> List[Dict]:
        """
        Generate dimer structures from donor + acceptors.

        CRITICAL: For 4-3 and 4-3-bridge crosslinks, the terminal D-Ala
        is cleaved from the donor during crosslinking.
        """
        dimers = []
        structure_column = columns["inferred"]["structure"]
        mass_column = columns["inferred"]["mass"]

        # Get lost residues for this crosslink type
        lost_residues = self._get_lost_residues(crosslink_type)

        # Calculate mass of lost residues.
        # Monoisotopic mass of free D-alanine (C3H7NO2), matching the
        # monoisotopic convention used by the "Theo (Da)" mass database.
        RESIDUE_MASSES = {"A": Decimal("89.047679")}  # D-Alanine
        lost_mass = sum(RESIDUE_MASSES.get(r, Decimal("0")) for r in lost_residues)

        # Water loss during peptide bond formation
        water_loss = Decimal("18.0106")

        # Remove oligomerization state from donor, then set aside any modification
        # tag (e.g. "gm-AEJAA (Anh)" -> stem "gm-AEJAA", tag " (Anh)") so trimming
        # below operates on the stem - modifications act on the glycan/side chains,
        # not the stem's terminal residue, and must be re-attached afterwards.
        donor_base = donor_structure.rsplit("|", 1)[0]
        donor_stem, donor_mod_tag = split_modification_tag(donor_base)

        # CRITICAL: Remove terminal D-Ala from donor when the crosslink type
        # loses one - the donor loses its terminal residue during transpeptidation.
        # This covers both pentapeptide->tetrapeptide trimming (4-3: "...AA" ->
        # "...A") and tetrapeptide->tripeptide trimming (3-3: "...A" -> "..."),
        # since both just remove exactly one trailing D-Ala.
        # Preserve bridges like [GGGGG] - don't trim immediately after a closing bracket.
        if "A" in lost_residues:
            if donor_stem.endswith("A") and not donor_stem.endswith("]A"):
                donor_stem = donor_stem[:-1]  # "gm-AQKAA" -> "gm-AQKA", "gm-AEJA" -> "gm-AEJ"

        donor_base = donor_stem + donor_mod_tag

        for _, acceptor_row in acceptors_df.iterrows():
            acceptor_structure = acceptor_row["structure"]
            acceptor_mass = Decimal(str(acceptor_row["mass"]))

            # Remove oligomerization state from acceptor
            acceptor_base = acceptor_structure.rsplit("|", 1)[0]

            # Calculate dimer mass (already accounts for lost_mass)
            dimer_mass = donor_mass + acceptor_mass - water_loss - lost_mass

            # Dimer structure notation: DONOR=ACCEPTOR[crosslink_type]|2
            # All crosslink types use = as separator; the bracket annotation
            # distinguishes the mechanism (4-3, 3-3, 4-3-bridge, Custom...).
            dimer_structure = f"{donor_base}={acceptor_base}[{crosslink_type}]|2"

            dimers.append({structure_column: dimer_structure, mass_column: float(dimer_mass)})

            LOGGER.debug(f"    Generated: {dimer_structure} (m/z: {float(dimer_mass):.4f})")

        LOGGER.info(f"  Generated {len(dimers)} dimers with {crosslink_type} crosslinks")

        return dimers


def build_dimers_for_species(
    matched_monomers_df: pd.DataFrame,
    theoretical_database_df: pd.DataFrame,
    species_code: Optional[str] = None,
    columns: dict = COLUMNS,
    strict_mode: bool = True,
    custom_rule: Optional[CustomCrosslinkRule] = None,
    donor_abundance_threshold: float = 0.9,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Convenience function to generate dimers using accurate biochemical rules.

    Parameters
    ----------
    matched_monomers_df : pd.DataFrame
        DataFrame of matched monomers from PGFinder (MS1 results)
    theoretical_database_df : pd.DataFrame
        DataFrame of all theoretical structures from database
    species_code : str, optional
        Species identifier. Exactly one of species_code/custom_rule must be given.
    columns : dict
        Column name mapping
    strict_mode : bool
        If True: Enforce exact literature rules (default)
        If False: Allow bridge variants for novel discovery
        Ignored when custom_rule is given.
    custom_rule : CustomCrosslinkRule, optional
        User-defined donor/acceptor regex rule. Exactly one of
        species_code/custom_rule must be given.
    donor_abundance_threshold : float
        Fraction (0.0-1.0) of the eligible donor pool's cumulative intensity
        to cover when selecting donors - see
        `AccurateDimerBuilder._find_eligible_donors_within_threshold`.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (theoretical dimers ready for matching, donors used to build them)
    """
    builder = AccurateDimerBuilder(
        species_code=species_code,
        custom_rule=custom_rule,
        strict_mode=strict_mode,
        donor_abundance_threshold=donor_abundance_threshold,
    )
    return builder.build_dimers_from_detected_monomers(matched_monomers_df, theoretical_database_df, columns)
