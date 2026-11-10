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
from plottingTrajectError import findXref, E
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

    def bananaPlot(self, cutoff=None):
        """Plot (R, x3) for each method and the reference in a 2x2 grid (tokamak only)."""
        names = list(self.results)                      # three methods, then "ref"
        fig, axs = plt.subplots(2, 2, figsize=(10, 10))

        for ax, name in zip(axs.ravel(), names):
            r = self.results[name]
            m = np.isfinite(r["X"]).all(axis=1)         # drop inf/nan (e.g. after a blow-up)
            if cutoff is not None:
                m &= r["t"] <= cutoff
            R = np.hypot(r["X"][m, 0], r["X"][m, 1])
            ax.plot(R, r["X"][m, 2], lw=0.5)
            ax.set_xlabel("$R$")
            ax.set_ylabel("$x_3$")
            ax.set_title("Reference" if name == "ref" else name)

        fig.suptitle(f"{self.sys['name']}: $(R, x_3)$ projection")
        fig.tight_layout()
        plt.show()
    def errorPlot(self, cutoff=None, stride=1):
        """Plot ||x_method(t_n) - x_ref(t_n)||_2 against t_n for all methods (log y-axis)."""
        ref = self.results["ref"]
        fig, ax = plt.subplots(figsize=(8, 5))

        for name, r in self.results.items():
            if name == "ref":
                continue
            err = np.linalg.norm(r["X"] - ref["X"], axis=1)
            err[~np.isfinite(err)] = np.nan              # blown-up points are skipped
            m = np.ones_like(err, dtype=bool)
            if cutoff is not None:
                m &= r["t"] <= cutoff
            ax.semilogy(r["t"][m][::stride], err[m][::stride], lw=0.8, label=name)

        ax.set_xlabel("$t_n$")
        ax.set_ylabel(r"$\|x^{(n)} - x_{\mathrm{ref}}(t_n)\|_2$")
        ax.set_title(f"{self.sys['name']}: error over time")
        ax.legend()
        fig.tight_layout()
        plt.show()

    def energyPlot(self, cutoff=None, stride=1):
        """Energy error |E_n - E_0| in a 2x2 grid: one panel per method, plus all methods together."""
        U = self.sys["U"]
        names = [n for n in self.results if n != "ref"]
        fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
        axs = axs.ravel()

        for k, name in enumerate(names):
            r = self.results[name]
            En = E(r["X"], r["V"], U)
            err = np.abs(En - En[0])
            err[~np.isfinite(err)] = np.nan              # skip blown-up points
            m = np.ones_like(err, dtype=bool)
            if cutoff is not None:
                m &= r["t"] <= cutoff
            t, e = r["t"][m][::stride], err[m][::stride]
            axs[k].semilogy(t, e, lw=0.6)                # own panel
            axs[3].semilogy(t, e, lw=1.2 if name == "RK4" else 0.6, label=name, zorder=3 if name == "RK4" else 1)   # combined panel
            axs[k].set_title(name)

        axs[3].set_title("All methods")
        
        axs[3].legend()
        for ax in axs:
            ax.set_xlabel("$t_n$")
        for ax in axs[::2]:
            ax.set_ylabel(r"$|E_n - E_0|$")

        fig.suptitle(f"{self.sys['name']}: energy error")
        fig.tight_layout()
        plt.show()

test = Traject("tokamak",T=1e3)
#test.trajectPlot()
test.energyPlot()
