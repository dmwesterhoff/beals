#!/usr/bin/env python3
"""Check that chunked certification results cover [2, Zimmert bound] exactly.

Usage: certify/check_results.py <field> [step ...]   (default steps: phase1 generators)

For each step, every results/<field>/<step>/*.out must end in "STATUS OK", use
the Zimmert bound recorded in results/<field>/setup.out, and the chunk ranges must be
contiguous from 2 to the bound. Exits non-zero on any gap, overlap or failure.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"


def setup_value(field_dir: Path, key: str) -> int:
    for line in (field_dir / "setup.out").read_text().splitlines():
        if line.startswith(key + " "):
            return int(line.split()[1])
    raise SystemExit(f"{key} missing from {field_dir}/setup.out")


def check_step(field_dir: Path, step: str, bound: int, prime_count: int) -> str:
    outs = sorted((field_dir / step).glob("*.out"))
    if not outs:
        raise SystemExit(f"{step}: no results")
    ranges: list[tuple[int, int]] = []
    ideals = 0
    primes = 0
    for path in outs:
        lines = path.read_text().splitlines()
        if not lines or lines[-1] != "STATUS OK":
            raise SystemExit(f"{step}: {path.name} incomplete")
        tag, a, b, bound_tag, bd = lines[0].split()
        if tag != "range" or bound_tag != "bound" or int(bd) != bound:
            raise SystemExit(f"{step}: {path.name} has wrong header {lines[0]!r}")
        ranges.append((int(a), int(b)))
        fields = lines[1].split()
        if fields[0] == "ideals":
            ideals += int(fields[1])
        elif fields[0] == "primes":
            primes += int(fields[1])
    ranges.sort()
    if ranges[0][0] != 2 or ranges[-1][1] != bound:
        raise SystemExit(f"{step}: covers [{ranges[0][0]}, {ranges[-1][1]}], expected [2, {bound}]")
    for (_, b1), (a2, _) in zip(ranges, ranges[1:], strict=False):
        if a2 != b1 + 1:
            raise SystemExit(f"{step}: gap or overlap between {b1} and {a2}")
    if primes and primes != prime_count:
        raise SystemExit(f"{step}: {primes} primes checked, expected pi(bound) = {prime_count}")
    extra = f", {ideals} prime ideals with verified generators" if ideals else f", {primes} primes"
    return f"{step.upper()} COMPLETE: {len(ranges)} contiguous chunks cover [2, {bound}]{extra}"


def main(argv: list[str]) -> int:
    if not argv:
        raise SystemExit(__doc__)
    field_dir = RESULTS / argv[0]
    steps = argv[1:] or ["phase1", "generators"]
    bound = setup_value(field_dir, "zimmert_bound")
    prime_count = setup_value(field_dir, "primepi_zimmert_bound")
    for step in steps:
        print(f"{argv[0]}: " + check_step(field_dir, step, bound, prime_count))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
