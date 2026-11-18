"""Render transparent-in-method, white-background SVG figures for the README."""

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
import os

import numpy as np

from filtered_variance import C2, C3_LOWER, EXTREMAL_HEAVY, GAMMA, paired_slice

BLUE = "#0072B2"
ORANGE = "#A74300"
TEAL = "#006B5B"
INK = "#172739"
MUTED = "#4D6173"
GRID = "#D6E0E7"


def identity_header(figure, title, subtitle):
    figure.text(0.035, 0.945, title, fontsize=17, fontweight="bold", color=INK, va="top")
    figure.text(0.035, 0.86, subtitle, fontsize=10.8, color=MUTED, va="top")


def interval(ax, a, b, y, color, label):
    ax.hlines(y, a, b, color=color, linewidth=3, zorder=2)
    ax.plot(a, y, marker="o", ms=7, mec=color, mew=1.7, mfc="white", zorder=3)
    ax.plot(b, y, marker="o", ms=7, mec=color, mew=1.2, mfc=color, zorder=3)
    ax.text((a+b)/2, y+0.18, label, ha="center", va="bottom", color=INK, fontsize=9.5)


def schedules(plt):
    fig, panels = plt.subplots(1, 2, figsize=(9.4, 4.0))
    fig.subplots_adjust(left=0.065, right=0.969, top=0.70, bottom=0.18, wspace=0.24)
    identity_header(fig, "Information intervals and reveal times",
                    r"The support is $(a,b]$. Open markers are issues and filled markers are maturities.")
    for ax in panels:
        ax.set_xlim(-0.42, 8.42)
        ax.set_ylim(0.22, 3.65)
        ax.set_xticks(range(0, 9, 2))
        ax.set_yticks([])
        ax.set_xlabel("Discrete time", labelpad=6)
        ax.spines[["left", "top", "right"]].set_visible(False)
    panels[0].set_title("Laminar: disjoint or nested", loc="left", fontsize=12, fontweight="bold")
    interval(panels[0], 0, 8, 3, BLUE, "(0, 8]")
    interval(panels[0], 0, 4, 2.20, BLUE, "(0, 4]")
    interval(panels[0], 4, 8, 1.91, TEAL, "(4, 8]")
    interval(panels[0], 0, 2, 1.22, BLUE, "(0, 2]")
    interval(panels[0], 2, 4, 0.93, TEAL, "(2, 4]")
    panels[0].axvline(2, ymin=0.10, ymax=0.31, color=MUTED, lw=1, ls="--")
    panels[0].text(2, 0.39, "maturity at next issue", fontsize=8.5, ha="center", color=MUTED)
    panels[1].set_title("Crossing: overlap without containment", loc="left",
                        fontsize=12, fontweight="bold")
    panels[1].axvspan(2, 5, ymin=0.20, ymax=0.83, color="#F9E9DE", zorder=0)
    interval(panels[1], 0, 5, 2.65, BLUE, "(0, 5]")
    interval(panels[1], 2, 8, 1.4, ORANGE, "(2, 8]")
    panels[1].text(3.5, 0.80, "common support (2, 5]", ha="center", fontsize=9, color=MUTED)
    return fig


def constants(plt):
    fig, ax = plt.subplots(figsize=(9.4, 3.65))
    fig.subplots_adjust(left=0.355, right=0.96, top=0.73, bottom=0.22)
    identity_header(fig, "Sharp factors and an unresolved constant",
                    r"Each row refers to a different schedule class. The bar marks proved bounds on $c_3$.")
    ax.set_ylim(-0.55, 2.46)
    ax.set_xlim(0.743, 0.887)
    ax.set_yticks([2, 1, 0], ["Laminar depth ≤ 2", "At most 3 two-leaf stars",
                               "Unrestricted laminar depth ≤ 3"])
    ax.set_xlabel(r"Worst-case ratio $V^\star/U$")
    ax.set_xticks(np.arange(0.75, 0.89, 0.025))
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.tick_params(length=0, axis="y", pad=9)
    ax.spines[["left", "top", "right"]].set_visible(False)
    ax.scatter(C2, 2, s=83, color=BLUE, zorder=2)
    ax.scatter(GAMMA, 1, s=82, marker="D", color=TEAL, zorder=2)
    ax.text(C2+0.003, 2, rf"$c_2\approx{C2:.4f}$", va="center", color=BLUE, fontsize=11)
    ax.text(GAMMA+0.0035, 1, rf"$\gamma\approx{GAMMA:.4f}$", va="center", color=TEAL, fontsize=11)
    ax.hlines(0, C3_LOWER, GAMMA, color=ORANGE, linewidth=4)
    ax.plot([C3_LOWER, GAMMA], [0, 0], linestyle="none", marker="|", ms=14, mew=1.8, color=ORANGE)
    ax.text(C3_LOWER, 0.17, r"$(1+4/\pi)/3$", color=ORANGE, ha="center", fontsize=10)
    ax.text(GAMMA, 0.17, r"$\gamma$", color=ORANGE, ha="center", fontsize=11)
    return fig


def weight_slice(plt):
    fig, ax = plt.subplots(figsize=(9.4, 4.45))
    fig.subplots_adjust(left=0.1, right=0.965, top=0.73, bottom=0.145)
    identity_header(fig, "A symmetric slice through the seven-weight optimum",
                    r"The same coefficients $(x,x,1,1)$ have different Gram norms under different pairings.")
    x = np.linspace(1, 6, 501)
    mixed = [paired_slice(value)[1] for value in x]
    like = [paired_slice(value, like_with_like=True)[1] for value in x]
    ax.plot(x, mixed, color=BLUE, lw=2.7, label="Heavy paired with unit")
    ax.plot(x, like, color=TEAL, lw=2.3, ls="--", label="Like paired with like")
    ax.axhline(GAMMA, color=ORANGE, lw=1.1, ls=":", zorder=1)
    ax.scatter([EXTREMAL_HEAVY], [GAMMA], s=70, facecolor="white", edgecolor=ORANGE, lw=2, zorder=4)
    ax.annotate(r"$x=(1+3\sqrt{2})/2$"
                "\n" r"$r=2,\quad V^\star/U=\gamma$",
                xy=(EXTREMAL_HEAVY, GAMMA), xytext=(3.48, 0.847),
                color=INK, va="center", fontsize=10,
                arrowprops={"arrowstyle": "->", "color": ORANGE})
    ax.set_xlim(1, 6)
    ax.set_ylim(0.806, 0.946)
    ax.set_xlabel(r"Heavy coefficient $x$ (lighter coefficient fixed at one)")
    ax.set_ylabel(r"$V^\star/U$ after minimizing over $r\geq0$")
    ax.legend(loc="upper left", frameon=False, fontsize=9.5)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("figures"))
    parser.add_argument("--preview", action="store_true", help="also write PNG inspection copies")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="fve-mpl-") as config:
        os.environ["MPLCONFIGDIR"] = config
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.rcParams.update({
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "font.size": 10.5,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "axes.edgecolor": MUTED,
            "xtick.color": MUTED,
            "ytick.color": INK,
            "axes.axisbelow": True,
            "svg.hashsalt": "filtered-variance-extremals-v1",
            "svg.fonttype": "path",
        })
        specs = [
            ("information-intervals", schedules, "Laminar and crossing deterministic information intervals."),
            ("sharp-factors", constants, "Two sharp factors and lower and upper bounds for unrestricted depth three."),
            ("seven-weight-slice", weight_slice, "Minimum variance quotient over the outer weight for two pairings."),
        ]
        for stem, make_plot, title in specs:
            figure = make_plot(plt)
            target = args.out_dir / f"{stem}.svg"
            figure.savefig(target, metadata={"Date": None, "Title": title})
            svg = target.read_text(encoding="utf-8")
            target.write_text("\n".join(line.rstrip() for line in svg.splitlines()) + "\n",
                              encoding="utf-8", newline="\n")
            if args.preview:
                figure.savefig(args.out_dir / f"{stem}.png", dpi=144)
            plt.close(figure)


if __name__ == "__main__":
    main()
