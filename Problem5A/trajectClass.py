import sys as _sys              # aliased, since `sys` is your system dict
from pathlib import Path

root = Path(__file__).resolve().parent.parent   # project root, one level above Problem5A
for p in (root, root / "symmetric", root / "Boris"):
    if str(p) not in _sys.path:
        _sys.path.insert(0, str(p))

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from create_systems import create_2d, create_tokamak
from RungeKutta4Method import RK4
from plottingTrajectError import findXref
from symmetricMethod import symmetricMultistep
from BorisMethod import boris




class Traject():
    def __init__(self,system,h=0.1,T=10e5):
        self.sys = create_2d() if system == "2d" else create_tokamak()
        self.h, self.T = h, T

        self.results = {}
        for method in (RK4, boris, symmetricMultistep):
            X, V, t, name = method(self.sys, h, T)
            print(f"{name} done")
            self.results[name] = dict(X=X, V=V, t=t)


        t = self.results["RK4"]["t"]
        self.results["ref"] = dict(X=findXref(self.sys, t), t=t)
        print(f"Ref Done")

    def trajectPlot(self, cutoff=None):
        """Plot (x1, x2) for each method in a 2x2 grid; the fourth panel is a 3D plot of the reference."""
        names = [n for n in self.results if n != "ref"]
        fig = plt.figure(figsize=(10, 10))

        for k, name in enumerate(names, start=1):
            r = self.results[name]
            m = np.isfinite(r["X"]).all(axis=1)          
            if cutoff is not None:
                m &= r["t"] <= cutoff
            ax = fig.add_subplot(2, 2, k)
            ax.plot(r["X"][m, 0], r["X"][m, 1], lw=0.5)
            ax.set_xlabel("$x_1$")
            ax.set_ylabel("$x_2$")
            ax.set_title(name)
            ax.set_aspect("equal")

        ref = self.results["ref"]
        m = np.isfinite(ref["X"]).all(axis=1)
        if cutoff is not None:
            m &= ref["t"] <= cutoff
        ax3 = fig.add_subplot(2, 2, 4, projection="3d")
        ax3.plot(ref["X"][m, 0], ref["X"][m, 1], ref["X"][m, 2], lw=0.5)
        ax3.set_xlabel("$x_1$")
        ax3.set_ylabel("$x_2$")
        ax3.set_zlabel("$x_3$")
        ax3.set_title("Reference (3D)")

        fig.suptitle(f"{self.sys['name']}: trajectories")
        fig.tight_layout()
        plt.show()

test = Traject("tokamak",T=1e3)
test.trajectPlot()
