#!/usr/bin/env python3
"""build_lab_meeting_deck.py ----
Assemble the self-contained lab-meeting HTML deck.

INPUT  : scripts/10_figures/lab_meeting_deck_template.html  (with __D01__..__D06__ placeholders)
         results/figures/12_deck/d0{1..6}_*.png
OUTPUT : results/figures/12_deck/lab_meeting_deck.html  (figures embedded as base64 data URIs)
Usage  : python3 scripts/10_figures/build_lab_meeting_deck.py

The output is a single file that opens in any browser with no external dependency
(offline-safe for projecting at the lab meeting). It is regenerated from the tracked
template + figures, so the large embedded file itself is gitignored.
"""
import base64
import pathlib
import sys

ROOT = pathlib.Path("/FAST/gr10634/gaozy/aml_niche_net")
TEMPLATE = ROOT / "scripts" / "10_figures" / "lab_meeting_deck_template.html"
FIGDIR = ROOT / "results" / "figures" / "12_deck"
OUT = FIGDIR / "lab_meeting_deck.html"

# placeholder -> figure file (stable, matches d_deck_figures.R outputs)
FIGURES = {
    "__D01__": "d01_cohort.png",
    "__D02__": "d02_label_quality.png",
    "__D03__": "d03_h1_null.png",
    "__D04__": "d04_h2.png",
    "__D05__": "d05_h3.png",
    "__D06__": "d06_ccc_content.png",
}


def data_uri(path: pathlib.Path) -> str:
    b64 = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def main() -> int:
    if not TEMPLATE.exists():
        print(f"ERROR: template not found: {TEMPLATE}", file=sys.stderr)
        return 1
    html = TEMPLATE.read_text(encoding="utf-8")

    for token, fname in FIGURES.items():
        fpath = FIGDIR / fname
        if not fpath.exists():
            print(f"ERROR: figure missing: {fpath}", file=sys.stderr)
            return 1
        html = html.replace(token, data_uri(fpath))

    # Fail loudly if any placeholder survived (typo / renamed figure)
    leftover = [t for t in FIGURES if t in html]
    if leftover:
        print(f"ERROR: unsubstituted placeholders: {leftover}", file=sys.stderr)
        return 1

    OUT.write_text(html, encoding="utf-8")
    kb = OUT.stat().st_size / 1024
    print(f"[deck] wrote {OUT}  ({kb:.0f} KB, {len(FIGURES)} figures embedded)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
