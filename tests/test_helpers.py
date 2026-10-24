"""Unit tests for the pure helper functions and the packet table.

These need no recordings or SDR hardware, so they always run.
"""
import re
from pathlib import Path

import numpy as np
import pytest

from helpers import (complex_mix, ecef_to_lla, fletcher_checksum,
                     lla_to_ecef, quad_interp, reverse_endian, rrcosfilter)
from orbcomm_packet import packet_dict

GOLD_DIR = Path(__file__).parent / "gold"


# --- checksum -------------------------------------------------------------

def _gold_full_packets():
    """Yield (hex_string, checksum_ok) for every full 12-byte packet printed
    in the gold files ('Unrecognized packet' lines show the whole packet)."""
    pattern = re.compile(r"^(### )?Unrecognized packet: ([0-9A-F]{24})$")
    for gold in sorted(GOLD_DIR.glob("*.txt")):
        for line in gold.read_text().splitlines():
            m = pattern.match(line.strip())
            if m:
                yield m.group(2), m.group(1) is None


def test_fletcher_checksum_matches_gold_packets():
    # A packet is valid when the checksum computed over the WHOLE packet
    # (including its own checksum bytes) is 0000.
    packets = list(_gold_full_packets())
    assert len(packets) > 50, "gold files should provide plenty of vectors"
    assert any(not ok for _, ok in packets), "need at least one failing vector"
    for hex_str, ok in packets:
        assert (fletcher_checksum(hex_str) == "0000") == ok, hex_str


def test_fletcher_checksum_simple_values():
    assert fletcher_checksum("00") == "0000"
    assert fletcher_checksum("01") == "0101"      # sum2=1, sum1=1
    assert fletcher_checksum("0101") == "0302"    # sum1=2, sum2=1+2


def test_fletcher_checksum_pads_odd_length():
    assert fletcher_checksum("ABC") == fletcher_checksum("ABC0")


# --- bit reversal ---------------------------------------------------------

@pytest.mark.parametrize("before, after", [
    ("01", "80"), ("80", "01"), ("F0", "0F"), ("1A", "58"), ("00", "00"),
])
def test_reverse_endian_known_bytes(before, after):
    assert reverse_endian(before) == after


def test_reverse_endian_is_its_own_inverse():
    for s in ("65A8F9", "1A5C0FFE", "000102030405"):
        assert reverse_endian(reverse_endian(s)) == s


# --- coordinate conversion -----------------------------------------------

def test_lla_to_ecef_known_points():
    x, y, z = lla_to_ecef(0.0, 0.0, 0.0)
    assert (x, y, z) == pytest.approx((6378137.0, 0.0, 0.0))
    x, y, z = lla_to_ecef(90.0, 0.0, 0.0)
    assert z == pytest.approx(6356752.314245, abs=1e-3)   # WGS84 semi-minor
    assert x == pytest.approx(0.0, abs=1e-3)


@pytest.mark.parametrize("lat", [-89, -45, 0, 30, 44.8, 89])
@pytest.mark.parametrize("lon", [-170, -116, 0, 90, 179])
@pytest.mark.parametrize("alt", [0, 1000, 715e3, 1200e3])
def test_ecef_lla_round_trip(lat, lon, alt):
    # ecef_to_lla is a closed-form approximation; measured worst case is
    # about 7e-8 deg and 9 mm, so these tolerances leave a safe margin.
    lat2, lon2, alt2 = ecef_to_lla(*lla_to_ecef(lat, lon, alt))
    assert lat2 == pytest.approx(lat, abs=1e-6)
    assert lon2 == pytest.approx(lon, abs=1e-6)
    assert alt2 == pytest.approx(alt, abs=0.05)


# --- DSP helpers ----------------------------------------------------------

def test_complex_mix_chunks_are_phase_coherent():
    rng = np.random.default_rng(0)
    x = (rng.standard_normal(1000) + 1j * rng.standard_normal(1000)).astype(np.complex64)
    fs, shift = 1.2288e6, 12345.0
    whole, _ = complex_mix(x, shift, fs)
    first, phase = complex_mix(x[:400], shift, fs)
    second, _ = complex_mix(x[400:], shift, fs, phase=phase)
    np.testing.assert_allclose(np.concatenate([first, second]), whole, atol=1e-5)


def test_complex_mix_shifts_a_tone():
    fs, n = 48000.0, 4800
    t = np.arange(n) / fs
    tone = np.exp(2j * np.pi * 1000.0 * t).astype(np.complex64)
    shifted, _ = complex_mix(tone, -1000.0, fs)
    # tone at +1 kHz mixed by -1 kHz lands at DC
    assert np.abs(np.mean(shifted)) == pytest.approx(1.0, abs=1e-3)


def test_quad_interp_returns_sub_sample_peak_offset():
    # parabola samples around a peak: offset from the max bin is 0.25
    assert quad_interp(np.array([0, 1, 4, 3.0, 0])) == pytest.approx(0.25)
    # symmetric neighbours -> peak exactly on the bin
    assert quad_interp(np.array([0, 3, 4, 3.0, 0])) == pytest.approx(0.0)


def test_rrcosfilter_shape():
    n = 101
    _, h = rrcosfilter(n, 0.35, 1 / 4800.0, 38400.0)
    assert len(h) == n
    assert np.all(np.isfinite(h))
    assert np.argmax(h) == n // 2
    np.testing.assert_allclose(h[1:], h[1:][::-1])   # symmetric about centre


# --- packet table ---------------------------------------------------------

HEX_PACKET_LEN = 24      # a standard packet is 12 bytes = 24 hex chars
EPHEMERIS_LEN = 48       # ephemeris packets are 24 bytes


def test_packet_headers_are_consistent():
    headers = set()
    for name, info in packet_dict.items():
        assert int(info["header"], 2) == int(info["hex_header"], 16), name
        assert info["hex_header"] not in headers, f"duplicate header: {name}"
        headers.add(info["hex_header"])


def test_packet_parts_fit_inside_packet():
    for name, info in packet_dict.items():
        limit = EPHEMERIS_LEN if name == "Ephemeris" else HEX_PACKET_LEN
        for part, (start, stop) in info["message_parts"]:
            assert 0 <= start < stop <= limit, f"{name}.{part}"
