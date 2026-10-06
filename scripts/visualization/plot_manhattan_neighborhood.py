#!/usr/bin/env python3
"""マンハッタン距離1～2の24格子点を正確に列挙して描く説明用の科学図。"""
import argparse
from itertools import product
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    points = [(x,y,z) for x,y,z in product(range(-2,3), repeat=3)
              if 1 <= abs(x)+abs(y)+abs(z) <= 2]
    assert len(points) == 24
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111, projection="3d")
    # 連続領域の境界は正八面体。探索対象はその中の離散的な格子点。
    faces = [[(2*sx,0,0),(0,2*sy,0),(0,0,2*sz)]
             for sx,sy,sz in product((-1,1), repeat=3)]
    ax.add_collection3d(Poly3DCollection(faces, alpha=0.08,
                                        facecolor="#4472C4", edgecolor="#4472C4"))
    for distance, color in ((1,"#23854B"),(2,"#4472C4")):
        subset=[p for p in points if sum(abs(v) for v in p)==distance]
        ax.scatter(*zip(*subset),s=55,c=color,depthshade=False,
                   label=f"Distance {distance}: {len(subset)} sites")
    ax.scatter(0,0,0,s=90,c="black",depthshade=False,label="Source endpoint")
    for setter in (ax.set_xlim,ax.set_ylim,ax.set_zlim): setter(-2.5,2.5)
    for setter in (ax.set_xticks,ax.set_yticks,ax.set_zticks): setter([-2,-1,0,1,2])
    ax.set_xlabel("dx"); ax.set_ylabel("dy"); ax.set_zlabel("dz")
    ax.set_box_aspect((1,1,1)); ax.view_init(elev=22,azim=-55)
    ax.set_title("Manhattan radius 2: 24 neighboring sites")
    ax.legend(loc="upper left",fontsize=9)
    fig.tight_layout(); args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=250,bbox_inches="tight",facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
