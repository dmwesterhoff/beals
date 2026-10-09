"""Tests for the L11 class group / unit certification scripts.

They drive the patched gp binary built by certify/build_pari.sh. Fast tests
only; the full Phase 1 / generator runs are done by certify/run_chunks.sh.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
GP = ROOT / ".tools" / "pari" / "bin" / "gp"

pytestmark = pytest.mark.skipif(not GP.exists(), reason="patched gp not built; run certify/build_pari.sh")


def gp(code: str, env: dict[str, str] | None = None, timeout: int = 600) -> str:
    """Run GP code (from stdin, one statement per line) in the repo root."""
    proc = subprocess.run(
        [str(GP), "-q", "-D", "parisizemax=1000000000"],
        input=code,
        capture_output=True,
        text=True,
        cwd=ROOT,
        env={**os.environ, **(env or {})},
        timeout=timeout,
    )
    return proc.stdout + proc.stderr


def gp_script(script: str, env: dict[str, str], timeout: int = 600) -> str:
    proc = subprocess.run(
        [str(GP), "-q", "-D", "parisizemax=1000000000", script],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        cwd=ROOT,
        env={**os.environ, **env},
        timeout=timeout,
    )
    return proc.stdout + proc.stderr


# Zimmert, Invent. Math. 62 (1981), Tabelle 1: (r1, r2, gamma, d^(1/2)/Na).
ZIMMERT_TABLE = [
    (2, 0, 2.41, 1.760),
    (1, 1, 1.84, 3.355),
    (1, 2, 1.20, 21.11),
    (0, 4, 0.87, 385.5),
    (10, 0, 0.56, 5.854e4),
    (20, 0, 0.36, 4.332e11),
    (0, 10, 0.46, 5.736e8),
    (0, 50, 0.18, 2.417e55),
]


@pytest.mark.parametrize("r1,r2,gam,expected", ZIMMERT_TABLE)
def test_zimmert_reproduces_published_table(r1, r2, gam, expected):
    out = gp(f'read("certify/zimmert.gp");\nprint(exp(zimmertZ({r1},{r2},{gam},zimmertalpha({gam}))));\n')
    value = float(out.strip().replace(" E", "e"))
    # Zimmert rounds his entries down to 4 significant digits.
    assert expected <= value < expected * 1.001


def test_zimmert_proof_form_matches_printed_satz2():
    out = gp(
        'read("certify/zimmert.gp");\n'
        "print(vecmax(vector(50, i, my(r1 = i % 7, r2 = i % 5 + 1, g = i / 37.);"
        " abs(zimmertZ(r1, r2, g, zimmertalpha(g))"
        " - zimmertZ_satz2(r1, r2, g, zimmertalpha(g))))) < 1e-30);\n"
    )
    assert out.strip() == "1"


def test_zimmert_bound_for_l11_is_pinned():
    out = gp(
        'read("certify/zimmert.gp");\n'
        "nf = nfinit(t^22 + 2*t^11 - 4); b = zimmertbest(2, 10);\n"
        "print(ceil(sqrt(abs(nf.disc)) * exp(-zimmertZ(2, 10, b[2], zimmertalpha(b[2])))) + 1);\n"
    )
    assert out.strip() == "152405135"


def test_hyperelliptic_algebra_is_l11():
    # C_{11,0}: 5Y^2 = (4 th + 12) X^11 - (4 th + 7) over K = Q(th), th^2 + th - 1 = 0.
    out = gp(
        "f = rnfequation(nfinit(y^2 + y - 1), (4*y + 12)*x^11 - (4*y + 7));\n"
        "print(nfisisom(f, t^22 + 2*t^11 - 4) != 0);\n"
    )
    assert out.strip() == "1"


@pytest.fixture
def class_number_two_bnf(tmp_path):
    path = tmp_path / "qsqrtm5.bnf"
    out = gp(f'B = bnfinit(t^2 + 5, 1); print(B.no); writebin("{path}", B);\n')
    assert out.strip() == "2"
    return path


def test_generators_rejects_nonprincipal_ideal(class_number_two_bnf, tmp_path):
    out = gp_script(
        "certify/generators.gp",
        {
            "BNF": str(class_number_two_bnf),
            "P1_A": "2",
            "P1_B": "50",
            "P1_BOUND": "100",
            "P1_OUT": str(tmp_path / "gens.out"),
        },
    )
    assert "nonprincipal" in out
    assert not (tmp_path / "gens.out").exists()


def test_generators_accepts_class_number_one_field(tmp_path):
    bnf = tmp_path / "q5.bnf"
    gp(f'writebin("{bnf}", bnfinit(t^2 - t - 1, 1));\n')
    out_file = tmp_path / "gens.out"
    gp_script(
        "certify/generators.gp",
        {"BNF": str(bnf), "P1_A": "2", "P1_B": "1000", "P1_BOUND": "1000", "P1_OUT": str(out_file)},
    )
    lines = out_file.read_text().splitlines()
    assert lines[-1] == "STATUS OK"


def test_charrank_detects_non_saturated_units():
    out = gp(
        'read("certify/units.gp");\n'
        "nf = nfinit(t^2 - 2); u = nfalgtobasis(nf, 1 + t);\n"
        "print(charrank(nf, [nfalgtobasis(nf, -1), u])[1]);\n"
        "print(charrank(nf, [nfalgtobasis(nf, -1), nfeltpow(nf, u, 2)], 1000, 60)[1]);\n",
        env={"UNITS_LIBRARY_ONLY": "1"},
    )
    assert out.split() == ["2", "1"]


def test_l11_units_are_two_saturated(tmp_path):
    out_file = tmp_path / "units.out"
    gp_script("certify/units.gp", {"UNITS_OUT": str(out_file)})
    text = out_file.read_text()
    assert "F2_rank 12" in text
    assert text.splitlines()[-1] == "STATUS OK"


def test_phase1_small_range_passes():
    out_file = ROOT / "results" / "test-phase1.tmp"
    try:
        gp_script(
            "certify/phase1.gp",
            {"P1_A": "100003", "P1_B": "100500", "P1_BOUND": "152405135", "P1_OUT": str(out_file)},
        )
        assert out_file.read_text().splitlines()[-1] == "STATUS OK"
    finally:
        out_file.unlink(missing_ok=True)
