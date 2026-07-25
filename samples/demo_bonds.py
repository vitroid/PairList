#!/usr/bin/env python3
"""Generate a demo figure of neighbor bonds (no GenIce required).

Usage:
    python samples/demo_bonds.py
    python samples/demo_bonds.py -o benchmark/demo.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

# Allow running from a checkout without installing.
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import pairlist as pl


def min_image_segments(
    pos_frac: np.ndarray, pairs: list[tuple[int, int, float]], cell: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return absolute (N_bond, 2, 2) xy segments and distances, using MIC."""
    segs = []
    dists = []
    for i, j, r in pairs:
        d = pos_frac[j] - pos_frac[i]
        d -= np.floor(d + 0.5)
        # skip bonds that wrap far in z-heavy systems for a clean xy view
        a = pos_frac[i] @ cell
        b = a + d @ cell
        segs.append([[a[0], a[1]], [b[0], b[1]]])
        dists.append(r)
    if not segs:
        return np.zeros((0, 2, 2)), np.zeros(0)
    return np.asarray(segs), np.asarray(dists)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("benchmark/demo.png"),
        help="output image path",
    )
    parser.add_argument("-n", type=int, default=120, help="number of particles")
    parser.add_argument("--box", type=float, default=12.0, help="cubic box length")
    parser.add_argument("--cutoff", type=float, default=2.2, help="pair cutoff")
    parser.add_argument("--seed", type=int, default=2)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    cell = np.eye(3) * args.box
    pos = rng.random((args.n, 3))
    xyz = pos @ cell

    pairs = list(pl.pairs_iter(pos, maxdist=args.cutoff, cell=cell))
    segs, dists = min_image_segments(pos, pairs, cell)

    # Style: slate points, teal→amber distance map on cool paper
    bg = "#eef2f5"
    fig, ax = plt.subplots(figsize=(6.2, 6.2), dpi=160)
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)

    if len(segs):
        lc = LineCollection(
            segs,
            array=dists,
            cmap="viridis",
            linewidths=0.95,
            alpha=0.88,
            zorder=1,
        )
        ax.add_collection(lc)
        cbar = fig.colorbar(lc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("pair distance", fontsize=9)
        cbar.ax.tick_params(labelsize=8)

    ax.scatter(
        xyz[:, 0],
        xyz[:, 1],
        s=30,
        c="#16324f",
        edgecolors=bg,
        linewidths=0.45,
        zorder=2,
    )

    # box outline (xy)
    L = args.box
    ax.plot([0, L, L, 0, 0], [0, 0, L, L, 0], color="#7a8a99", lw=1.0, zorder=0)

    ax.set_aspect("equal")
    ax.set_xlim(-0.4, L + 0.4)
    ax.set_ylim(-0.4, L + 0.4)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        f"pairlist demo — N={args.n}, cutoff={args.cutoff} (xy, PBC)",
        fontsize=11,
        pad=10,
    )
    for spine in ax.spines.values():
        spine.set_color("#b7c2cc")
    ax.tick_params(colors="#4a5560", labelsize=8)
    ax.xaxis.label.set_color("#4a5560")
    ax.yaxis.label.set_color("#4a5560")
    ax.title.set_color("#1b2430")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {args.output} ({len(pairs)} pairs)")


if __name__ == "__main__":
    main()
