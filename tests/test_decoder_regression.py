"""Regression test: decode each recording and compare with the gold output.

For every tests/gold/<name>.txt this runs `file_decoder.py data/<name>.mat`
as a subprocess and compares the parsed result, not the raw text:

  * packet lines (type, hex payload, checksum-failed flag): exact match
  * header values (bit offset, packet counts, ...): exact match
  * residual frequency offset: within +/- FREQ_OFFSET_TOL_HZ
  * the satellite's own reported position: tight tolerance
  * values computed from our TLE propagator (the 'Ephem ...' lines and the
    differences): looser tolerance, so the pyephem -> skyfield swap can be
    checked against the same gold files. TUNE THESE when switching.

The .mat recordings are stored with git-LFS. If a recording is missing or is
an un-fetched LFS pointer, its test is skipped (run `git lfs pull`).
"""
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
GOLD_DIR = Path(__file__).parent / "gold"
DATA_DIR = REPO / "data"

FREQ_OFFSET_TOL_HZ = 2.0
# absolute tolerances for propagator-derived numbers
EPHEM_TOL = {"lat": 0.05, "lon": 0.05, "alt_km": 1.0, "vel_ms": 5.0,
             "dist_km": 2.0, "dvel_ms": 5.0}

HEADER_EXACT_KEYS = [
    "Satellites in recording", "SDR Sample rate", "SDR Center frequency",
    "Satellite frequencies", "Number of samples", "Number of possible packets",
    "Bit stream offset", "Number of packets decoded",
]
FREQ_RE = re.compile(r"Remaining frequency offset after doppler compensation: (-?[\d.]+) Hz")
SUMMARY_RE = re.compile(r"(\d+) packets with errors, (\d+) packets corrected, PER:\s+([\d.]+)%")
NUM = r"(-?\d+\.\d+)"
SAT_POS_RE = re.compile(rf"Sat Lat/Lon:\s+{NUM}, {NUM}, Altitude:\s+{NUM} km, Velocity:\s+{NUM} m/s")
EPH_POS_RE = re.compile(rf"Ephem Lat/Lon:\s+{NUM}, {NUM}, Altitude:\s+{NUM} km, Velocity:\s+{NUM} m/s")
DIST_RE = re.compile(rf"Difference in reported and ephemeris position:\s+{NUM} km")
DVEL_RE = re.compile(rf"Difference in reported and ephemeris velocity:\s+{NUM} m/s")


def parse_output(text):
    """Split decoder stdout into the pieces we compare."""
    lines = text.splitlines()
    header = {}
    for line in lines:
        if line.startswith("List of packets"):
            break
        key, _, value = line.partition(":")
        header[key.strip()] = value.strip()

    m = FREQ_RE.search(text)
    freq_offset = float(m.group(1)) if m else None

    packets, details = [], []      # details: indented lines belonging to a packet
    in_list = False
    for line in lines:
        if line.startswith("List of packets"):
            in_list = True
            continue
        if not in_list or SUMMARY_RE.match(line):
            continue
        if line.startswith("\t"):
            details.append(line.strip())
        elif line.strip():
            packets.append(line.rstrip())

    summary = SUMMARY_RE.search(text)
    return {
        "header": {k: header.get(k) for k in HEADER_EXACT_KEYS},
        "freq_offset": freq_offset,
        "packets": packets,
        "details": details,
        "summary": summary.groups() if summary else None,
    }


def _first(regex, details):
    for d in details:
        m = regex.search(d)
        if m:
            return [float(g) for g in m.groups()]
    return None


def _data_file_or_skip(name):
    mat = DATA_DIR / f"{name}.mat"
    if not mat.exists():
        pytest.skip(f"{mat} not found")
    with open(mat, "rb") as f:
        if f.read(40).startswith(b"version https://git-lfs"):
            pytest.skip(f"{mat} is an un-fetched git-LFS pointer (run: git lfs pull)")
    return mat


@pytest.fixture(scope="module")
def decode(tmp_path_factory):
    """Run file_decoder.py once per recording and cache the parsed result."""
    cache = {}

    def run(name):
        if name not in cache:
            mat = _data_file_or_skip(name)
            workdir = tmp_path_factory.mktemp(name)   # packets.txt lands here
            env = dict(os.environ, MPLBACKEND="Agg")
            proc = subprocess.run(
                [sys.executable, str(REPO / "file_decoder.py"), str(mat)],
                cwd=workdir, env=env, capture_output=True, text=True, timeout=600)
            assert proc.returncode == 0, proc.stderr[-2000:]
            cache[name] = parse_output(proc.stdout)
        return cache[name]

    return run


GOLD_NAMES = sorted(p.stem for p in GOLD_DIR.glob("*.txt"))


def test_gold_files_present():
    assert GOLD_NAMES, "no gold files found in tests/gold"


@pytest.mark.parametrize("name", GOLD_NAMES)
def test_header_values(decode, name):
    got = decode(name)
    gold = parse_output((GOLD_DIR / f"{name}.txt").read_text())
    assert got["header"] == gold["header"]


@pytest.mark.parametrize("name", GOLD_NAMES)
def test_frequency_offset(decode, name):
    got = decode(name)
    gold = parse_output((GOLD_DIR / f"{name}.txt").read_text())
    assert got["freq_offset"] == pytest.approx(gold["freq_offset"], abs=FREQ_OFFSET_TOL_HZ)


@pytest.mark.parametrize("name", GOLD_NAMES)
def test_packets_match(decode, name):
    got = decode(name)
    gold = parse_output((GOLD_DIR / f"{name}.txt").read_text())
    assert len(got["packets"]) == len(gold["packets"])
    assert got["packets"] == gold["packets"]
    assert got["summary"] == gold["summary"]


@pytest.mark.parametrize("name", GOLD_NAMES)
def test_ephemeris_details(decode, name):
    got = decode(name)
    gold = parse_output((GOLD_DIR / f"{name}.txt").read_text())
    assert len(got["details"]) == len(gold["details"])

    # The satellite's own reported time and position come from the packet
    # data, not from our propagator, so they must match tightly.
    gold_time = [d for d in gold["details"] if d.startswith("Current satellite time")]
    got_time = [d for d in got["details"] if d.startswith("Current satellite time")]
    assert got_time == gold_time
    g, w = _first(SAT_POS_RE, got["details"]), _first(SAT_POS_RE, gold["details"])
    assert (g is None) == (w is None)
    if w:
        assert g == pytest.approx(w, abs=0.01)

    # Propagator-derived values: looser tolerances (see EPHEM_TOL).
    for regex, keys in [(EPH_POS_RE, ["lat", "lon", "alt_km", "vel_ms"]),
                        (DIST_RE, ["dist_km"]), (DVEL_RE, ["dvel_ms"])]:
        g, w = _first(regex, got["details"]), _first(regex, gold["details"])
        assert (g is None) == (w is None)
        if w:
            for gv, wv, key in zip(g, w, keys):
                assert gv == pytest.approx(wv, abs=EPHEM_TOL[key]), key
