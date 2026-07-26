#!/usr/bin/env python3
"""fig_cooption_composition.py ----
Stacked-bar figure: sender-compartment COMPOSITION of the top conserved CCC pathways,
AML vs healthy. Visualization only -- no conclusion text is baked into the figure
(interpretation + the within-study caveat live in CANDIDATE_M10_cooption_lens.md).

Reads the cohort-mean composition produced by cooption_redistribution.py and collapses
the 7 sender bins into three compartments:
  primitive = HSC_MPP, LMPP_GMP        (LSC-like / stem-progenitor)
  immune    = Mono_DC, T_NK, B_Plasma
  other     = Erythroid, Megakaryocyte

INPUT  : results/tables/05_ccc/cooption/sender_composition.csv
OUTPUT : results/figures/13_cooption/cooption_sender_composition.png
Usage  : python3 scripts/10_figures/fig_cooption_composition.py

Palette: dataviz validated categorical slots 1/2/3 (blue/orange/aqua) -- CVD-safe as a set.
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

ROOT = pathlib.Path("/FAST/gr10634/gaozy/aml_niche_net")
INP = ROOT / "results" / "tables" / "05_ccc" / "cooption" / "sender_composition.csv"
OUTDIR = ROOT / "results" / "figures" / "13_cooption"
OUT = OUTDIR / "cooption_sender_composition.png"

# validated palette (dataviz palette.md, light surface)
PRIM_C, IMM_C, OTH_C = "#eb6834", "#2a78d6", "#1baf7a"   # slot 2 / 1 / 3
SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID = "#e4e3df"

MACRO = {
    "HSC_MPP": "primitive", "LMPP_GMP": "primitive",
    "Mono_DC": "immune", "T_NK": "immune", "B_Plasma": "immune",
    "Erythroid": "other", "Megakaryocyte": "other",
}
# pathway display order: the two clear movers first, MIF (near-flat control) last
PATH_ORDER = ["LGALS9", "CD99", "MIF"]
PATH_LABEL = {"LGALS9": "Galectin-9\n(LGALS9)", "CD99": "CD99", "MIF": "MIF"}
GROUPS = ["AML", "healthy"]
N_LABEL = {"AML": "AML\n(n=121)", "healthy": "healthy\n(n=27)"}


def macro_shares() -> pd.DataFrame:
    c = pd.read_csv(INP)
    c["macro"] = c["sender_bin"].map(MACRO)
    m = (c.groupby(["ligand", "group", "macro"])["mean_sender_share"].sum()
         .reset_index())
    piv = m.pivot_table(index=["ligand", "group"], columns="macro",
                        values="mean_sender_share").fillna(0.0)
    for col in ("primitive", "immune", "other"):
        if col not in piv.columns:
            piv[col] = 0.0
    return piv


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    piv = macro_shares()

    fig, ax = plt.subplots(figsize=(10.5, 6.2), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    bw = 0.34
    offs = {"AML": -0.5 * bw - 0.03, "healthy": 0.5 * bw + 0.03}
    xticks, xticklabels = [], []

    for i, lig in enumerate(PATH_ORDER):
        for grp in GROUPS:
            x = i + offs[grp]
            row = piv.loc[(lig, grp)]
            p, im, ot = row["primitive"] * 100, row["immune"] * 100, row["other"] * 100
            # stack: primitive (baseline-anchored; it is the quantity of interest) -> immune -> other
            ax.bar(x, p, bw, bottom=0, color=PRIM_C, edgecolor=SURFACE, linewidth=2, zorder=3)
            ax.bar(x, im, bw, bottom=p, color=IMM_C, edgecolor=SURFACE, linewidth=2, zorder=3)
            ax.bar(x, ot, bw, bottom=p + im, color=OTH_C, edgecolor=SURFACE, linewidth=2, zorder=3)
            # direct-label the primitive share only (the story quantity), white centered
            if p >= 8:
                ax.text(x, p / 2, f"{p:.0f}%", ha="center", va="center",
                        color="white", fontsize=12.5, fontweight="bold", zorder=4)
            ax.text(x, -4.5, N_LABEL[grp], ha="center", va="top",
                    color=INK2, fontsize=10.5, linespacing=1.05)
        xticks.append(i)
        xticklabels.append(PATH_LABEL[lig])

    # pathway labels sit below the AML/healthy labels
    ax.set_xticks(xticks)
    ax.set_xticklabels([""] * len(xticks))
    for i, lig in enumerate(PATH_ORDER):
        ax.text(i, -14.5, PATH_LABEL[lig], ha="center", va="top",
                color=INK, fontsize=13.5, fontweight="bold", linespacing=1.0)

    ax.set_ylim(0, 100)
    ax.set_xlim(-0.62, len(PATH_ORDER) - 0.38)
    ax.set_ylabel("Share of pathway signal by sender (%)", fontsize=13, color=INK)
    ax.set_yticks(range(0, 101, 20))
    ax.tick_params(axis="y", labelsize=11.5, colors=INK2, length=0)
    ax.tick_params(axis="x", length=0)
    ax.yaxis.grid(True, color=GRID, linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(GRID)

    ax.set_title("Who sends the top conserved pathways: sender composition, AML vs healthy",
                 fontsize=15.5, color=INK, fontweight="bold", pad=26, loc="left", x=0.0)
    ax.text(0.0, 1.045, "Cohort mean share of each pathway's significant signal, by sender compartment.",
            transform=ax.transAxes, fontsize=10.5, color=INK2, ha="left")

    # short labels only (compartment membership is defined in the doc caption); placed in the
    # right margin so it never overlaps the bars, which reach 100%.
    legend = [
        Patch(facecolor=PRIM_C, label="Primitive"),
        Patch(facecolor=IMM_C, label="Immune"),
        Patch(facecolor=OTH_C, label="Other"),
    ]
    ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(1.015, 1.0), frameon=False,
              fontsize=11.5, handlelength=1.2, handleheight=1.2, labelcolor=INK,
              borderaxespad=0.0, labelspacing=0.7)

    # room for the two-tier x labels (bottom) and the right-margin legend
    fig.subplots_adjust(bottom=0.17, top=0.86, left=0.085, right=0.86)
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"[fig] wrote {OUT}  ({OUT.stat().st_size/1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
