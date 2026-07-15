"""
mzML reader for MS2 spectra.

Extracts MS2 scans from standard mzML files (Thermo, Bruker, Agilent, Waters).
Requires pyteomics: install with  pip install pgfinder[ms2]

Usage
-----
    from pgfinder.ms2.io import iter_ms2_scans, read_ms2_scans

    # Memory-efficient — yields one scan at a time
    for scan in iter_ms2_scans("experiment.mzML"):
        print(scan["precursor_mz"], scan["peaks"].shape)

    # Load everything into a list
    scans = read_ms2_scans("experiment.mzML")
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import pandas as pd


def iter_ms2_scans(path: str | Path) -> Iterator[dict]:
    """
    Yield MS2 scan dicts from an mzML file, one at a time.

    Each yielded dict contains:

        scan_id          str           raw scan identifier from the file
        rt               float         retention time in minutes
        precursor_mz     float         selected-ion m/z of the precursor
        precursor_charge int           precursor charge state (0 if unknown)
        peaks            pd.DataFrame  columns: "mz" (float64), "intensity" (float64)

    MS1 scans and any MS2 scan that has no precursor or no peaks are silently
    skipped.  The caller gets only clean, usable MS2 spectra.

    Parameters
    ----------
    path : str or Path
        Path to the .mzML file.

    Raises
    ------
    ImportError
        If pyteomics is not installed (``pip install pgfinder[ms2]``).
    FileNotFoundError
        If the file does not exist.
    """
    try:
        from pyteomics import mzml
    except ImportError as exc:
        raise ImportError(
            "pyteomics is required to read mzML files. " "Install it with:  pip install pgfinder[ms2]"
        ) from exc

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"mzML file not found: {path}")

    with mzml.MzML(str(path)) as reader:
        for spectrum in reader:
            if spectrum.get("ms level") != 2:
                continue

            mz_arr = spectrum.get("m/z array")
            int_arr = spectrum.get("intensity array")
            if mz_arr is None or int_arr is None or len(mz_arr) == 0:
                continue

            # Precursor
            prec_list = spectrum.get("precursorList", {})
            precursors = prec_list.get("precursor", [])
            if not precursors:
                continue

            selected_ions = precursors[0].get("selectedIonList", {}).get("selectedIon", [{}])
            ion = selected_ions[0] if selected_ions else {}

            precursor_mz = ion.get("selected ion m/z")
            if precursor_mz is None:
                continue

            charge = ion.get("charge state", 0)

            # Retention time — pyteomics returns a unitfloat; plain float() strips units
            scan_entries = spectrum.get("scanList", {}).get("scan", [{}])
            rt = float(scan_entries[0].get("scan start time", float("nan"))) if scan_entries else float("nan")

            yield {
                "scan_id": spectrum.get("id", ""),
                "rt": rt,
                "precursor_mz": float(precursor_mz),
                "precursor_charge": int(charge) if charge else 0,
                "peaks": pd.DataFrame({"mz": mz_arr, "intensity": int_arr}),
            }


def read_ms2_scans(path: str | Path) -> list[dict]:
    """
    Read all MS2 scans from an mzML file into a list.

    For large files (>1 GB) prefer ``iter_ms2_scans`` to avoid loading
    every scan into memory at once.

    Parameters
    ----------
    path : str or Path
        Path to the .mzML file.

    Returns
    -------
    list[dict]
        List of scan dicts as described in ``iter_ms2_scans``.
    """
    return list(iter_ms2_scans(path))
