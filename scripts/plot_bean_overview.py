"""Plot one mid-pressure equilibrium from each indentation family."""
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as Polygon
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from eq_process import read_geqdsk


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    args = parser.parse_args()
    root = args.dataset
    rows = json.loads((root/"catalog.json").read_text())
    shapes = sorted({r["indent_id"] for r in rows})
    fig, axes = plt.subplots(2, 4, figsize=(13, 7.8), layout="constrained")
    for ax, shape in zip(axes.flat, shapes):
        family = [r for r in rows if r["indent_id"] == shape]
        row = family[len(family)//2]
        d = read_geqdsk(root/row["gfile"])
        rr, zz = np.meshgrid(d.rgrid, d.zgrid, indexing="ij")
        inside = Polygon(np.column_stack((d.rbbbs, d.zbbbs))).contains_points(
            np.column_stack((rr.ravel(), zz.ravel()))).reshape(rr.shape)
        u = np.ma.masked_where(~inside, (d.psirz-d.simag)/(d.sibry-d.simag))
        ax.contour(rr, zz, u, levels=[.05, .2, .4, .6, .8, .95], cmap="viridis", linewidths=1)
        ax.plot(np.r_[d.rbbbs, d.rbbbs[0]], np.r_[d.zbbbs, d.zbbbs[0]], color="#222222", lw=1.4)
        ax.plot(d.rmaxis, d.zmaxis, "+", color="#b91c1c", ms=7)
        ax.set_title(f"{shape}: i = {row['indentation']:g}\n"
                     f"p0 = {row['p0_kPa']:g} kPa, beta_av = {row['beta_av_percent']:.2f}%", fontsize=10)
        ax.set(xlabel="R [m]", ylabel="Z [m]", aspect="equal")
        ax.grid(alpha=.15)
    fig.suptitle("Bean equilibrium archive: LCFS, normalized poloidal-flux contours and magnetic axis", fontsize=13)
    fig.savefig(root/"overview.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
