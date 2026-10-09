"""Tests for the free-software reimplementation of Stoll's criterion (stoll/criterion.gp)."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
GP = ROOT / ".tools" / "pari" / "bin" / "gp"

pytestmark = pytest.mark.skipif(not GP.exists(), reason="patched gp not built; run certify/build_pari.sh")


def gp_fn(body: str, timeout: int = 600) -> list[str]:
    """Define and call a GP function `tst()` with the given body; return printed lines."""
    script = 'read("stoll/criterion.gp");\ntst() =\n{\n' + body + "\n}\ntst();\nquit;\n"
    path = ROOT / "results" / "stoll" / "_test.gp"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(script)
    try:
        proc = subprocess.run(
            [str(GP), "-q", "-D", "parisizemax=1000000000", str(path)],
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            cwd=ROOT,
            timeout=timeout,
        )
    finally:
        path.unlink(missing_ok=True)
    errors = [ln for ln in proc.stderr.splitlines() if "***" in ln and "Warning" not in ln]
    assert not errors, proc.stderr
    return proc.stdout.split()


def test_l2sqrt_squares_and_nonsquares():
    out = gp_fn(
        "  my(l = 7, lam = Mod(t, t^7 - 2), a = 3 + lam + lam^5, s);\n"
        "  s = l2sqrt(a^2, l, 80); print(lval(s^2 - a^2, l) >= 70 * l);\n"
        "  print(l2sqrt(Mod(5, t^7 - 2), l, 80) == 0);\n"  # unramified class: not a square
        "  print(l2sqrt(lam, l, 80) == 0);\n"  # odd valuation
        "  print(l2sqrt(Mod(1 + 8 * 3, t^7 - 2), l, 80) != 0);"  # 25 is a square
    )
    assert out == ["1", "1", "1", "1"]


@pytest.mark.parametrize("ell", [7, 11, 13])
def test_square_class_basis_has_full_dimension(ell):
    out = gp_fn(f"  my(nf = nfinit(t^{ell} - 2), P = idealprimedec(nf, 2)[1]); print(#sqclassbasis(nf, P));")
    assert out == [str(ell + 2)]


@pytest.mark.parametrize("ell", [7, 11, 13])
def test_halving_matches_stoll_conjecture_8_6(ell):
    out = gp_fn(
        f"  my(l = {ell}, nf = nfinit(t^l - 2), P = idealprimedec(nf, 2)[1], B = sqclassbasis(nf, P));\n"
        "  my(v = [sigmaclass(l, 5, 120)[1], sigmaclass(l, -3, 120)[1], sigmaconj(l)]);\n"
        "  v = apply(z -> hilbvec(nf, lift(z), B, P), v);\n"
        '  print(v[1] == v[2], " ", v[1] == v[3], " ", v[1] != vector(#B));'
    )
    assert out == ["1", "1", "1"]


def test_criterion_passes_for_l7_and_detects_classes_in_image():
    out = gp_fn(
        "  my(R = stollcriterion(7, bnfinit(t^7 - 2, 1)), th = Mod(-t^5 / 2, t^7 - 2));\n"
        '  print(mapget(R, "PASS"));\n'
        # mu((1,1)) is a global Selmer element, so its class at 2 lies in r(S): must be flagged
        "  print(inimage2(R, (-4/5) * (th - 1)));\n"
        '  print(inimage2(R, mapget(R, "nf").pol * 0 + 1));'  # the trivial class is in the image
    )
    assert out == ["1", "1", "1"]


def test_criterion_rejects_wieferich_and_small_l():
    out = gp_fn("  print(iferr(stollcriterion(5, bnfinit(t^5 - 2, 1)); 0, E, 1));")
    assert out == ["1"]
