"""
Tests for pgfinder.ms2.io — mzML reader.

The fixture mzML is generated inline (no external file needed): three spectra
containing one MS1 scan (which must be skipped) and two MS2 scans.  Binary
arrays are zlib-compressed 64-bit little-endian floats, base64-encoded — the
standard mzML encoding.
"""

from __future__ import annotations

import textwrap

import pytest

from pgfinder.ms2.io import iter_ms2_scans, read_ms2_scans

# ---------------------------------------------------------------------------
# Fixture: minimal valid mzML with 1 MS1 + 2 MS2 scans
# ---------------------------------------------------------------------------

# Pre-computed base64(zlib(pack('<Nd', *values))) for the test arrays.
# Reproduce with:
#   import base64, struct, zlib
#   def b64(fs): return base64.b64encode(zlib.compress(struct.pack(f'<{len(fs)}d',*fs))).decode()
_MS1_MZ_B64 = "eJxjYACBSAcwxZAJoQ8UOQAAFFgCtQ=="  # [100.0, 200.0, 300.0]
_MS1_INT_B64 = "eJxjYAACh34HBjA9H0rXOwAAIFoDLg=="  # [1000.0, 2000.0, 500.0]
_S1_MZ_B64 = "eJzrdEx4eoEpyOFN4A651uYwh/v+vdPznNIdBCIst5xoynQAAO+gDhY="  # [72.044, 90.055, 186.076, 204.087]
_S1_INT_B64 = "eJxjYAACh34HBjBdD6XnQ+iCdgcAO2gEZQ=="  # [1000.0, 500.0, 2000.0, 750.0]
_S2_MZ_B64 = "eJwTiLDccqIp0+H2z7qsPZ/KHABBaQhp"  # [204.087, 367.171]
_S2_INT_B64 = "eJxjYACCjs0OIIqhYLkDABIEAtM="  # [5000.0, 3000.0]

_MINIMAL_MZML = textwrap.dedent(
    f"""\
    <?xml version="1.0" encoding="utf-8"?>
    <mzML xmlns="http://psi.hupo.org/ms/mzml"
          xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
          xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd">
      <cvList count="2">
        <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology"
            version="4.1.30" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>
        <cv id="UO" fullName="Unit Ontology" version="09:04:2014"
            URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>
      </cvList>
      <fileDescription>
        <fileContent>
          <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
        </fileContent>
      </fileDescription>
      <softwareList count="1">
        <software id="sw" version="1.0">
          <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value=""/>
        </software>
      </softwareList>
      <dataProcessingList count="1">
        <dataProcessing id="dp">
          <processingMethod order="0" softwareRef="sw">
            <cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" value=""/>
          </processingMethod>
        </dataProcessing>
      </dataProcessingList>
      <run>
        <spectrumList count="3" defaultDataProcessingRef="dp">

          <!-- MS1 scan — must be skipped by iter_ms2_scans -->
          <spectrum index="0" id="scan=1" defaultArrayLength="3">
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
            <scanList count="1">
              <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
              <scan>
                <cvParam cvRef="MS" accession="MS:1000016" name="scan start time"
                         value="0.5" unitAccession="UO:0000031" unitName="minute" unitCvRef="UO"/>
              </scan>
            </scanList>
            <binaryDataArrayList count="2">
              <binaryDataArray encodedLength="32">
                <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>
                <cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>
                <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
                <binary>{_MS1_MZ_B64}</binary>
              </binaryDataArray>
              <binaryDataArray encodedLength="32">
                <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>
                <cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>
                <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
                <binary>{_MS1_INT_B64}</binary>
              </binaryDataArray>
            </binaryDataArrayList>
          </spectrum>

          <!-- MS2 scan 1: precursor 471.207 z=2, rt=1.0 min -->
          <spectrum index="1" id="scan=2" defaultArrayLength="4">
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>
            <scanList count="1">
              <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
              <scan>
                <cvParam cvRef="MS" accession="MS:1000016" name="scan start time"
                         value="1.0" unitAccession="UO:0000031" unitName="minute" unitCvRef="UO"/>
              </scan>
            </scanList>
            <precursorList count="1">
              <precursor>
                <selectedIonList count="1">
                  <selectedIon>
                    <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z"
                             value="471.207" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2"/>
                  </selectedIon>
                </selectedIonList>
              </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
              <binaryDataArray encodedLength="56">
                <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>
                <cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>
                <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
                <binary>{_S1_MZ_B64}</binary>
              </binaryDataArray>
              <binaryDataArray encodedLength="44">
                <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>
                <cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>
                <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
                <binary>{_S1_INT_B64}</binary>
              </binaryDataArray>
            </binaryDataArrayList>
          </spectrum>

          <!-- MS2 scan 2: precursor 942.415 z=1, rt=2.5 min -->
          <spectrum index="2" id="scan=3" defaultArrayLength="2">
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>
            <scanList count="1">
              <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
              <scan>
                <cvParam cvRef="MS" accession="MS:1000016" name="scan start time"
                         value="2.5" unitAccession="UO:0000031" unitName="minute" unitCvRef="UO"/>
              </scan>
            </scanList>
            <precursorList count="1">
              <precursor>
                <selectedIonList count="1">
                  <selectedIon>
                    <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z"
                             value="942.415" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="1"/>
                  </selectedIon>
                </selectedIonList>
              </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
              <binaryDataArray encodedLength="32">
                <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>
                <cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>
                <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
                <binary>{_S2_MZ_B64}</binary>
              </binaryDataArray>
              <binaryDataArray encodedLength="28">
                <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>
                <cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>
                <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
                <binary>{_S2_INT_B64}</binary>
              </binaryDataArray>
            </binaryDataArrayList>
          </spectrum>

        </spectrumList>
      </run>
    </mzML>
"""
)


@pytest.fixture
def mzml_path(tmp_path):
    """Write the minimal test mzML to a temp file and return its Path."""
    p = tmp_path / "test.mzML"
    p.write_text(_MINIMAL_MZML, encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


def test_file_not_found():
    """iter_ms2_scans must raise FileNotFoundError for a missing file."""
    with pytest.raises(FileNotFoundError):
        list(iter_ms2_scans("/does/not/exist.mzML"))


# ---------------------------------------------------------------------------
# MS1 filtering
# ---------------------------------------------------------------------------


def test_ms1_scans_are_skipped(mzml_path):
    """The fixture has 1 MS1 scan — it must not appear in the output."""
    scans = read_ms2_scans(mzml_path)
    assert len(scans) == 2


# ---------------------------------------------------------------------------
# Returned dict keys and types
# ---------------------------------------------------------------------------


def test_scan_keys(mzml_path):
    """Each scan dict must contain the expected keys."""
    scans = read_ms2_scans(mzml_path)
    required = {"scan_id", "rt", "precursor_mz", "precursor_charge", "peaks"}
    for scan in scans:
        assert required.issubset(set(scan.keys()))


def test_peaks_dataframe_columns(mzml_path):
    """peaks must be a DataFrame with exactly 'mz' and 'intensity' columns."""
    import pandas as pd

    scans = read_ms2_scans(mzml_path)
    for scan in scans:
        assert isinstance(scan["peaks"], pd.DataFrame)
        assert set(scan["peaks"].columns) == {"mz", "intensity"}


def test_peaks_nonempty(mzml_path):
    """Every returned scan must have at least one peak."""
    scans = read_ms2_scans(mzml_path)
    for scan in scans:
        assert len(scan["peaks"]) > 0


# ---------------------------------------------------------------------------
# Precursor values
# ---------------------------------------------------------------------------


def test_precursor_mz_scan1(mzml_path):
    """First MS2 scan has precursor m/z = 471.207."""
    scans = read_ms2_scans(mzml_path)
    assert scans[0]["precursor_mz"] == pytest.approx(471.207, abs=1e-3)


def test_precursor_charge_scan1(mzml_path):
    """First MS2 scan has charge state 2."""
    scans = read_ms2_scans(mzml_path)
    assert scans[0]["precursor_charge"] == 2


def test_precursor_mz_scan2(mzml_path):
    """Second MS2 scan has precursor m/z = 942.415."""
    scans = read_ms2_scans(mzml_path)
    assert scans[1]["precursor_mz"] == pytest.approx(942.415, abs=1e-3)


def test_precursor_charge_scan2(mzml_path):
    """Second MS2 scan has charge state 1."""
    scans = read_ms2_scans(mzml_path)
    assert scans[1]["precursor_charge"] == 1


# ---------------------------------------------------------------------------
# Retention time
# ---------------------------------------------------------------------------


def test_retention_time_scan1(mzml_path):
    """First MS2 scan rt = 1.0 min."""
    scans = read_ms2_scans(mzml_path)
    assert scans[0]["rt"] == pytest.approx(1.0, abs=1e-4)


def test_retention_time_scan2(mzml_path):
    """Second MS2 scan rt = 2.5 min."""
    scans = read_ms2_scans(mzml_path)
    assert scans[1]["rt"] == pytest.approx(2.5, abs=1e-4)


# ---------------------------------------------------------------------------
# Peak m/z values decoded correctly
# ---------------------------------------------------------------------------


def test_peaks_mz_scan1(mzml_path):
    """First MS2 scan peaks include the expected m/z values."""
    scans = read_ms2_scans(mzml_path)
    mz_vals = scans[0]["peaks"]["mz"].tolist()
    assert mz_vals == pytest.approx([72.044, 90.055, 186.076, 204.087], abs=1e-3)


def test_peaks_mz_scan2(mzml_path):
    """Second MS2 scan peaks include the expected m/z values."""
    scans = read_ms2_scans(mzml_path)
    mz_vals = scans[1]["peaks"]["mz"].tolist()
    assert mz_vals == pytest.approx([204.087, 367.171], abs=1e-3)


# ---------------------------------------------------------------------------
# iter_ms2_scans (generator) vs read_ms2_scans (list)
# ---------------------------------------------------------------------------


def test_iter_gives_same_results_as_read(mzml_path):
    """iter_ms2_scans and read_ms2_scans must produce identical output."""
    from_iter = list(iter_ms2_scans(mzml_path))
    from_read = read_ms2_scans(mzml_path)
    assert len(from_iter) == len(from_read)
    for a, b in zip(from_iter, from_read):
        assert a["precursor_mz"] == pytest.approx(b["precursor_mz"])
        assert a["precursor_charge"] == b["precursor_charge"]
        assert a["rt"] == pytest.approx(b["rt"])
