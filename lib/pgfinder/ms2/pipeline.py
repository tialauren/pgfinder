"""
MS2 scoring pipeline: match mzML scans against a PG mass library.

For each MS2 scan:
  1. Compute the neutral precursor mass from (m/z, charge).
  2. Find library structures whose theoretical mass matches within ppm tolerance.
  3. Predict HCD fragments for each candidate and score against the observed spectrum.
  4. Return all matches above min_score, one row per (scan, structure) pair.

Usage
-----
    from pgfinder.ms2.pipeline import score_mzml_against_library
    import pandas as pd

    library = pd.read_csv("e_coli_monomers_complex.csv")
    results = score_mzml_against_library("experiment.mzML", library)
    print(results.head())
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

PROTON = 1.007276  # Da


def _neutral_mass(precursor_mz: float, charge: int) -> float:
    return precursor_mz * charge - charge * PROTON


def _candidates(
    neutral_mass: float,
    library: pd.DataFrame,
    mass_col: str,
    ppm: float,
) -> pd.DataFrame:
    delta = library[mass_col] - neutral_mass
    ppm_error = (delta / neutral_mass).abs() * 1e6
    return library[ppm_error <= ppm].copy().assign(_ppm_error=ppm_error[ppm_error <= ppm])


def score_mzml_against_library(
    mzml_path: str | Path,
    library: pd.DataFrame,
    structure_col: str = "Structure",
    mass_col: str = "Monoisotopic Mass",
    precursor_tolerance_ppm: float = 10.0,
    fragment_tolerance_da: float = 0.02,
    charge_max: int = 3,
    reduced: bool = True,
    min_score: float = 0.0,
    charges_to_try: tuple[int, ...] = (1, 2, 3),
) -> pd.DataFrame:
    """
    Score every MS2 scan in an mzML file against a PG mass library.

    Parameters
    ----------
    mzml_path : str or Path
        Path to the .mzML file.
    library : pd.DataFrame
        Mass library with at least `structure_col` and `mass_col` columns.
        Accepts pgfinder's built-in format ("Structure", "Monoisotopic Mass")
        or ftrs output ("Inferred structure", "Theo (Da)").
    structure_col : str
        Name of the column containing structure strings.
    mass_col : str
        Name of the column containing monoisotopic neutral masses (Da).
    precursor_tolerance_ppm : float
        Mass tolerance for precursor matching (default 10 ppm).
    fragment_tolerance_da : float
        Mass tolerance for fragment matching in the scorer (default 0.02 Da).
    charge_max : int
        Maximum fragment charge state to predict (default 3).
    reduced : bool
        True (default) for NaBH4-reduced spectra; False for non-reduced (Kwan 2024).
    min_score : float
        Minimum cosine similarity score to include in output (default 0.0).
    charges_to_try : tuple[int]
        Charge states to try when a scan reports charge 0 (unknown).

    Returns
    -------
    pd.DataFrame
        One row per (scan, matched structure). Columns:
        scan_id, rt, precursor_mz, precursor_charge, structure,
        score, ppm_error.
        Sorted by score descending.
    """
    from pgfinder.ms2.io import iter_ms2_scans
    from pgfinder.ms2.scorer import score_spectrum

    rows = []

    for scan in iter_ms2_scans(mzml_path):
        charge = scan["precursor_charge"]
        mz = scan["precursor_mz"]
        peaks = scan["peaks"]

        # If charge is unknown, try several
        charges = [charge] if charge > 0 else list(charges_to_try)

        for z in charges:
            nm = _neutral_mass(mz, z)
            candidates = _candidates(nm, library, mass_col, precursor_tolerance_ppm)
            if candidates.empty:
                continue

            for _, row in candidates.iterrows():
                structure = row[structure_col]
                ppm_err = row["_ppm_error"]

                result = score_spectrum(
                    structure=structure,
                    experimental_peaks=peaks,
                    min_da_tolerance=fragment_tolerance_da,
                    charge_max=charge_max,
                    reduced=reduced,
                )
                if result is None or result["score"] < min_score:
                    continue

                rows.append(
                    {
                        "scan_id": scan["scan_id"],
                        "rt": scan["rt"],
                        "precursor_mz": mz,
                        "precursor_charge": z,
                        "structure": structure,
                        "score": result["score"],
                        "matched_peaks": result["matched_peaks"],
                        "total_theo": result["total_theo"],
                        "coverage": result["coverage"],
                        "ppm_error": round(ppm_err, 3),
                    }
                )

    if not rows:
        return pd.DataFrame(
            columns=["scan_id", "rt", "precursor_mz", "precursor_charge", "structure", "score", "matched_peaks", "total_theo", "coverage", "ppm_error"]
        )

    return pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
