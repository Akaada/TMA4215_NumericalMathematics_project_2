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
    """Class for solving a system with our methods. Stored in one object so that it is easier to work with in a longer notebook"""
    def __init__(self, system, h=0.1, T=1e5, only_symmetric=False):
        """Init of a traject class that solves the system with the necessary methods"""
        self.sys = create_2d() if system == "2d" else create_tokamak() # Create the system

        # Store necessary values 
        self.h, self.T = h, T
        self.only_symmetric = only_symmetric
        self.results = {}

        methods = (symmetricMultistep,) if only_symmetric else (RK4, boris, symmetricMultistep) # define our methods
        for method in methods: # Go trough the methods
            X, V, t, name = method(self.sys, h, T)
            print(f"{name} done")
            self.results[name] = dict(X=X, V=V, t=t) # Store results in a dictionary

        if not only_symmetric: # This is mainly so that we can test symmetric multistep over a longer time period without needing long runtime
            t = self.results["RK4"]["t"]
            self.results["ref"] = dict(X=findXref(self.sys, t), t=t) # Find a reference solution with an order 8 RungeKutta method
            print("Ref done")

    def trajectPlot(self, cutoff=None):
        """Plot (x1, x2) for each method in a 2x2 grid; the fourth panel is a 3D plot of the reference."""
        names = [n for n in self.results if n != "ref"] # Get the names that is not the reference solution
        fig = plt.figure(figsize=(10, 10)) # Init our figure

        for k, name in enumerate(names, start=1): # Loop over the solutions created by the methods and plot the trajectory
            r = self.results[name]
            m = np.isfinite(r["X"]).all(axis=1) # Make it so that if there is a blow up we still get a plot and we can see how it diverges         
            if cutoff is not None:  # So that we can plot for shorter time period
                m &= r["t"] <= cutoff
            ax = fig.add_subplot(2, 2, k)
            ax.plot(r["X"][m, 0], r["X"][m, 1], lw=0.5)
            ax.set_xlabel("$x_1$")
            ax.set_ylabel("$x_2$")
            ax.set_title(name)
            ax.set_aspect("equal")

        ref = self.results["ref"]
        m = np.isfinite(ref["X"]).all(axis=1) # Make it so that if there is a blow up we still get a plot and we can see how it diverges   
        if cutoff is not None: # So that we can plot for shorter time period
            m &= ref["t"] <= cutoff
        ax3 = fig.add_subplot(2, 2, 4, projection="3d") # Add the reference solution as a 3D-plot
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
            if cutoff is not None: # So that we can plot for a shorter time period
                m &= r["t"] <= cutoff 
            R = np.hypot(r["X"][m, 0], r["X"][m, 1]) # Calculate the Hypotenuse
            ax.plot(R, r["X"][m, 2], lw=0.5)
            ax.set_xlabel("$R$")
            ax.set_ylabel("$x_3$")
            ax.set_title("Reference" if name == "ref" else name) # set title, we store refrence as ref, so we change it for the visuals

        fig.suptitle(f"{self.sys['name']}: $(R, x_3)$ projection")
        fig.tight_layout()
        plt.show()
    def errorPlot(self, cutoff=None, stride=1):
        """Plot ||x_method(t_n) - x_ref(t_n)||_2 against t_n for all methods (log y-axis)."""
        ref = self.results["ref"] # Get the reference
        fig, ax = plt.subplots(figsize=(8, 5))

        for name, r in self.results.items(): # loop over results
            if name == "ref": # do not need reference
                continue
            err = np.linalg.norm(r["X"] - ref["X"], axis=1)
            err[~np.isfinite(err) | (err > 1e10)] = np.nan   # treat blow-up as missing              
            m = np.ones_like(err, dtype=bool)
            if cutoff is not None: # so we can plot for shorter time period
                m &= r["t"] <= cutoff
            ax.semilogy(r["t"][m][::stride], err[m][::stride], lw=0.8, label=name) # plot all methods in one plot

        ax.set_xlabel("$t_n$")
        ax.set_ylabel(r"$\|x^{(n)} - x_{\mathrm{ref}}(t_n)\|_2$")
        ax.set_title(f"{self.sys['name']}: error over time")
        ax.legend()
        fig.tight_layout()
        plt.show()

    def _energyError(self, name, cutoff, stride): # create a new one to work better inside the object
        """Return (t, |E_n - E_0|) for one method, with blow-up masked."""
        r = self.results[name]
        En = E(r["X"], r["V"], self.sys["U"]) # Calculate the energy
        err = np.abs(En - En[0]) # calculate the energy error
        err[~np.isfinite(err) | (err > 1e10)] = np.nan          # skip blown-up points
        m = np.ones_like(err, dtype=bool)
        if cutoff is not None: # so we can plot for shorter time period
            m &= r["t"] <= cutoff
        return r["t"][m][::stride], err[m][::stride]
    
    def energyPlot(self, cutoff=None, stride=1):
        """Energy error |E_n - E_0|: one panel if only_symmetric, otherwise a 2x2 grid."""
        if self.only_symmetric: # If we are only looking at symmetric we need a slightly different plot, so we just do it inside this if statement.
            name = "Symmetric Multistep"
            t, e = self._energyError(name, cutoff, stride)
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.semilogy(t, e, lw=0.6)
            ax.set_xlabel("$t_n$")
            ax.set_ylabel(r"$|E_n - E_0|$")
            ax.set_title(f"{self.sys['name']}: energy error, {name}, $h = {self.h}$")
            fig.tight_layout()
            plt.show()
            return

        names = [n for n in self.results if n != "ref"]
        fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
        axs = axs.ravel()

        for k, name in enumerate(names): # Loop over all the methods and plot once in its own panel, and once in the joint panel
            t, e = self._energyError(name, cutoff, stride)
            axs[k].semilogy(t, e, lw=0.6)                                    # own panel
            axs[3].semilogy(t, e, lw=1.2 if name == "RK4" else 0.6, label=name,
                            zorder=3 if name == "RK4" else 1)                # combined panel, need a higher zorder for RK4 since it usually is less turbulent, and will be lost if not
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